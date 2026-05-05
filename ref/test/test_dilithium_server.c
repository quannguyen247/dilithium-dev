#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <limits.h>

/* Platform-specific socket headers - must come before ../sign.h to avoid macro conflicts */
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <ws2tcpip.h>
  #include <windows.h>
  #define close(sock) closesocket(sock)
  #define ssize_t int
  /* Undefine potential macro conflicts that Windows headers define */
  #undef N
  #undef D
  #undef L
#else
  #include <sys/types.h>
  #include <sys/socket.h>
  #include <netinet/in.h>
  #include <arpa/inet.h>
  #include <unistd.h>
  #include <signal.h>
#endif

/* Dilithium headers - included after socket headers to avoid macro conflicts */
#include "../randombytes.h"
#include "../sign.h"

/* Configuration */
#define SERVER_PORT 5000
#define CHALLENGE_MAX 8192

/* Forward declarations */
static uint64_t get_time_ms(void);
static int send_all(int sock, const uint8_t *buf, size_t len);
static int recv_all(int sock, uint8_t *buf, size_t len);
static int send_blob(int sock, const uint8_t *data, uint32_t data_len);
static int recv_blob(int sock, uint8_t *buf, uint32_t buf_size, uint32_t *out_len);
static int receive_public_key(int sock, uint8_t *pk);
static int send_challenge(int sock, const uint8_t *challenge, size_t challenge_len);
static int receive_signature(int sock, uint8_t *signature, size_t *sig_len);
static int handle_client(int client_sock);

/* Get current time in milliseconds */
static uint64_t get_time_ms(void) {
#ifdef _WIN32
  return (uint64_t)GetTickCount64();
#else
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
    return 0;
  }
  return (uint64_t)ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
#endif
}

static int send_all(int sock, const uint8_t *buf, size_t len) {
  size_t total = 0;
  while (total < len) {
    ssize_t sent = send(sock, (const char *)buf + total, (int)(len - total), 0);
    if (sent < 0) {
#ifndef _WIN32
      if (errno == EINTR) {
        continue;
      }
#endif
      return -1;
    }
    if (sent == 0) {
      errno = ECONNRESET;
      return -1;
    }
    total += (size_t)sent;
  }
  return 0;
}

static int recv_all(int sock, uint8_t *buf, size_t len) {
  size_t total = 0;
  while (total < len) {
    ssize_t recvd = recv(sock, (char *)buf + total, (int)(len - total), 0);
    if (recvd == 0) {
      errno = ECONNRESET;
      return -1;
    }
    if (recvd < 0) {
#ifndef _WIN32
      if (errno == EINTR) {
        continue;
      }
#endif
      return -1;
    }
    total += (size_t)recvd;
  }
  return 0;
}

static int send_blob(int sock, const uint8_t *data, uint32_t data_len) {
  uint32_t len_net = htonl(data_len);
  if (send_all(sock, (const uint8_t *)&len_net, sizeof(len_net)) < 0) {
    return -1;
  }
  if (data_len == 0) {
    return 0;
  }
  return send_all(sock, data, data_len);
}

static int recv_blob(int sock, uint8_t *buf, uint32_t buf_size, uint32_t *out_len) {
  uint32_t len_net = 0;
  if (recv_all(sock, (uint8_t *)&len_net, sizeof(len_net)) < 0) {
    return -1;
  }

  uint32_t len = ntohl(len_net);
  if (len > buf_size) {
    errno = EMSGSIZE;
    return -1;
  }

  if (len > 0 && recv_all(sock, buf, len) < 0) {
    return -1;
  }

  *out_len = len;
  return 0;
}


/* Receive public key from client */
static int receive_public_key(int sock, uint8_t *pk) {
  printf("[*] Waiting for public key from client...\n");

  uint32_t pk_size = 0;
  if (recv_blob(sock, pk, CRYPTO_PUBLICKEYBYTES, &pk_size) < 0) {
    perror("recv() public key failed");
    return -1;
  }

  if (pk_size != (uint32_t)CRYPTO_PUBLICKEYBYTES) {
    fprintf(stderr, "Invalid public key size: %u (expected %d)\n", pk_size, CRYPTO_PUBLICKEYBYTES);
    return -1;
  }

  printf("[+] Public key received successfully (size: %u bytes)\n\n", pk_size);
  return 0;
}

/* Send challenge message to client */
static int send_challenge(int sock, const uint8_t *challenge, size_t challenge_len) {
  printf("[*] Sending challenge to client (size: %zu bytes)...\n", challenge_len);

  if (challenge_len > UINT32_MAX) {
    fprintf(stderr, "Challenge too large: %zu bytes\n", challenge_len);
    return -1;
  }

  if (send_blob(sock, challenge, (uint32_t)challenge_len) < 0) {
    perror("send() challenge failed");
    return -1;
  }

  printf("[+] Challenge sent successfully\n\n");
  return 0;
}

/* Receive signature from client */
static int receive_signature(int sock, uint8_t *signature, size_t *sig_len) {
  printf("[*] Waiting for signature from client...\n");

  uint32_t size = 0;
  if (recv_blob(sock, signature, CRYPTO_BYTES, &size) < 0) {
    perror("recv() signature failed");
    return -1;
  }

  *sig_len = (size_t)size;
  printf("[+] Signature received successfully (size: %zu bytes)\n\n", *sig_len);
  return 0;
}


static int handle_client(int client_sock) {
  uint64_t total_start = get_time_ms();

  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t challenge[CHALLENGE_MAX];
  size_t challenge_len = 0;
  uint8_t signature[CRYPTO_BYTES];
  size_t sig_len = 0;

  /* ============ STAGE 1: Receive Public Key from Client ============ */
  printf("[STAGE 1] Receiving public key from client...\n");

  uint64_t recv_pk_start = get_time_ms();
  if (receive_public_key(client_sock, pk) < 0) {
    fprintf(stderr, "Failed to receive public key\n");
    return 1;
  }
  uint64_t recv_pk_end = get_time_ms();

  printf("- Client's Public Key Size: %d bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("[+] Receive PK time: %llu ms\n\n",
         (unsigned long long)(recv_pk_end - recv_pk_start));

  /* ============ STAGE 2: Load Challenge & Send ============ */
  printf("[STAGE 2] Loading challenge and sending to client...\n");

  FILE *fin = fopen("test/input.txt", "rb");
  if (!fin) {
    fin = fopen("input.txt", "rb");
  }

  if (!fin) {
    const char *default_msg = "This is a test challenge message";
    size_t default_len = strlen(default_msg);
    memcpy(challenge, default_msg, default_len);
    challenge_len = default_len;
    printf("[WARNING] Cannot open input file, using default challenge\n");
  } else {
    challenge_len = fread(challenge, 1, CHALLENGE_MAX, fin);
    fclose(fin);
    if (challenge_len == 0) {
      const char *default_msg = "This is a test challenge message";
      size_t default_len = strlen(default_msg);
      memcpy(challenge, default_msg, default_len);
      challenge_len = default_len;
      printf("[WARNING] Empty input file, using default challenge\n");
    }
  }

  printf("- Challenge loaded: %zu bytes\n", challenge_len);

  uint64_t send_challenge_start = get_time_ms();
  if (send_challenge(client_sock, challenge, challenge_len) < 0) {
    fprintf(stderr, "Failed to send challenge\n");
    return 1;
  }
  uint64_t send_challenge_end = get_time_ms();

  printf("[+] Send challenge time: %llu ms\n\n",
         (unsigned long long)(send_challenge_end - send_challenge_start));

  /* ============ STAGE 3: Receive Signature ============ */
  printf("[STAGE 3] Receiving signature from client...\n");

  uint64_t recv_sig_start = get_time_ms();
  if (receive_signature(client_sock, signature, &sig_len) < 0) {
    fprintf(stderr, "Failed to receive signature\n");
    return 1;
  }
  uint64_t recv_sig_end = get_time_ms();

  printf("- Signature size: %zu bytes\n", sig_len);
  printf("[+] Receive signature time: %llu ms\n\n",
         (unsigned long long)(recv_sig_end - recv_sig_start));

  /* ============ STAGE 4: Verify Signature ============ */
  printf("[STAGE 4] Verifying signature...\n");

  uint64_t verify_start = get_time_ms();
  int verify_result = crypto_sign_verify(signature, sig_len, challenge, challenge_len,
                                         NULL, 0, pk);
  uint64_t verify_end = get_time_ms();

  printf("- Verification result: %s\n", verify_result == 0 ? "VALID" : "INVALID");
  printf("[+] Verification time: %llu ms\n\n",
         (unsigned long long)(verify_end - verify_start));

  /* ============ TIMING SUMMARY ============ */
  uint64_t total_end = get_time_ms();

  printf("===================================\n");
  printf("[TIMING SUMMARY]\n");
  printf("===================================\n");
  printf("Receive PK Time:           %llu ms\n",
         (unsigned long long)(recv_pk_end - recv_pk_start));
  printf("Send Challenge Time:       %llu ms\n",
         (unsigned long long)(send_challenge_end - send_challenge_start));
  printf("Receive Signature Time:    %llu ms\n",
         (unsigned long long)(recv_sig_end - recv_sig_start));
  printf("Verification Time:         %llu ms\n",
         (unsigned long long)(verify_end - verify_start));
  printf("-----------------------------------\n");
  printf("Total Time (from start):   %llu ms\n",
         (unsigned long long)(total_end - total_start));
  printf("===================================\n\n");

  printf("[KEY INFORMATION]\n");
  printf("- Client Public Key Size:  %d bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("- Signature Size:          %zu bytes\n", sig_len);
  printf("- Challenge Size:          %zu bytes\n", challenge_len);
  printf("===================================\n\n");

  printf("[+] Signature verification %s.\n", verify_result == 0 ? "OK" : "FAILED");

  return verify_result == 0 ? 0 : 1;
}

int main(void) {
  int listen_sock = -1;

  printf("\n========== Dilithium Server ==========\n");
  printf("Listening on port %d\n", SERVER_PORT);
  printf("======================================\n\n");

  /* Windows socket initialization */
#ifdef _WIN32
  WSADATA wsa_data;
  if (WSAStartup(MAKEWORD(2, 2), &wsa_data) != 0) {
    fprintf(stderr, "WSAStartup failed\n");
    return 1;
  }
#else
  signal(SIGCHLD, SIG_IGN);
  signal(SIGPIPE, SIG_IGN);
#endif

  /* ============ STAGE 0: Create Socket & Listen ============ */
  listen_sock = socket(AF_INET, SOCK_STREAM, 0);
  if (listen_sock < 0) {
    perror("socket() failed");
#ifdef _WIN32
    WSACleanup();
#endif
    return 1;
  }

  /* Allow socket address reuse */
  int reuse = 1;
  if (setsockopt(listen_sock, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(reuse)) < 0) {
    perror("setsockopt() failed");
    close(listen_sock);
#ifdef _WIN32
    WSACleanup();
#endif
    return 1;
  }

  struct sockaddr_in server_addr;
  memset(&server_addr, 0, sizeof(server_addr));
  server_addr.sin_family = AF_INET;
  server_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  server_addr.sin_port = htons(SERVER_PORT);

  if (bind(listen_sock, (struct sockaddr *)&server_addr, sizeof(server_addr)) < 0) {
    perror("bind() failed");
    close(listen_sock);
#ifdef _WIN32
    WSACleanup();
#endif
    return 1;
  }

  if (listen(listen_sock, SOMAXCONN) < 0) {
    perror("listen() failed");
    close(listen_sock);
#ifdef _WIN32
    WSACleanup();
#endif
    return 1;
  }

  printf("[*] Waiting for client connections...\n");

  while (1) {
    struct sockaddr_in client_addr;
    socklen_t client_addr_len = sizeof(client_addr);
    int client_sock = accept(listen_sock, (struct sockaddr *)&client_addr, &client_addr_len);
    if (client_sock < 0) {
#ifndef _WIN32
      if (errno == EINTR) {
        continue;
      }
#endif
      perror("accept() failed");
      continue;
    }

    printf("[+] Client connected from %s:%d\n\n", inet_ntoa(client_addr.sin_addr),
           ntohs(client_addr.sin_port));

#ifndef _WIN32
    pid_t pid = fork();
    if (pid == 0) {
      close(listen_sock);
      handle_client(client_sock);
      close(client_sock);
      _exit(0);
    }

    if (pid < 0) {
      perror("fork() failed");
      close(client_sock);
      continue;
    }

    close(client_sock);
#else
    handle_client(client_sock);
    close(client_sock);
#endif
  }

  close(listen_sock);
#ifdef _WIN32
  WSACleanup();
#endif
  return 0;
}


