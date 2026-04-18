#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Platform-specific socket headers - must come before ../sign.h to avoid macro conflicts */
#ifdef _WIN32
  #include <winsock2.h>
  #include <ws2tcpip.h>
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
#endif

/* Dilithium headers - included after socket headers to avoid macro conflicts */
#include "../randombytes.h"
#include "../sign.h"

/* Configuration */
#define SERVER_PORT 5000
#define BUFFER_SIZE 8192
#define MLEN 1200  /* Challenge size */

/* Forward declarations */
uint64_t get_time_ms(void);
int receive_public_key(int sock, uint8_t *pk);
int send_challenge(int sock, const uint8_t *challenge, size_t challenge_len);
int receive_signature(int sock, uint8_t *signature, size_t *sig_len);

/* Get current time in milliseconds */
uint64_t get_time_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
}


/* Receive public key from client */
int receive_public_key(int sock, uint8_t *pk) {
  printf("[*] Waiting for public key from client...\n");
  
  uint32_t pk_size;
  ssize_t recvd = recv(sock, (char *)&pk_size, sizeof(pk_size), 0);
  if (recvd < 0) {
    perror("recv() public key size failed");
    return -1;
  }
  
  if (recvd == 0) {
    fprintf(stderr, "Connection closed by client\n");
    return -1;
  }
  
  if (pk_size != (uint32_t)CRYPTO_PUBLICKEYBYTES) {
    fprintf(stderr, "Invalid public key size: %u (expected %d)\n", pk_size, CRYPTO_PUBLICKEYBYTES);
    return -1;
  }
  
  recvd = recv(sock, (char *)pk, CRYPTO_PUBLICKEYBYTES, 0);
  if (recvd < 0) {
    perror("recv() public key failed");
    return -1;
  }
  
  printf("[+] Public key received successfully (size: %zu bytes)\n\n", (size_t)recvd);
  return 0;
}

/* Send challenge message to client */
int send_challenge(int sock, const uint8_t *challenge, size_t challenge_len) {
  printf("[*] Sending challenge to client (size: %zu bytes)...\n", challenge_len);
  
  uint32_t size = (uint32_t)challenge_len;
  if (send(sock, (char *)&size, sizeof(size), 0) < 0) {
    perror("send() challenge size failed");
    return -1;
  }
  
  ssize_t sent = send(sock, (char *)challenge, challenge_len, 0);
  if (sent < 0) {
    perror("send() challenge failed");
    return -1;
  }
  
  printf("[+] Challenge sent successfully\n\n");
  return 0;
}

/* Receive signature from client */
int receive_signature(int sock, uint8_t *signature, size_t *sig_len) {
  printf("[*] Waiting for signature from client...\n");
  
  uint32_t size;
  ssize_t recvd = recv(sock, (char *)&size, sizeof(size), 0);
  if (recvd < 0) {
    perror("recv() signature size failed");
    return -1;
  }
  
  if (recvd == 0) {
    fprintf(stderr, "Connection closed by client\n");
    return -1;
  }
  
  if (size > BUFFER_SIZE) {
    fprintf(stderr, "Signature size too large: %u\n", size);
    return -1;
  }
  
  recvd = recv(sock, (char *)signature, size, 0);
  if (recvd < 0) {
    perror("recv() signature failed");
    return -1;
  }
  
  *sig_len = (size_t)recvd;
  printf("[+] Signature received successfully (size: %zu bytes)\n\n", *sig_len);
  return 0;
}


int main(void) {
  uint64_t total_start = get_time_ms();
  
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];  /* Client's public key */
  uint8_t challenge[BUFFER_SIZE];
  size_t challenge_len;
  uint8_t signature[BUFFER_SIZE];
  size_t sig_len;
  
  int listen_sock = -1;
  int client_sock = -1;
  
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
  
  if (listen(listen_sock, 1) < 0) {
    perror("listen() failed");
    close(listen_sock);
    #ifdef _WIN32
      WSACleanup();
    #endif
    return 1;
  }
  
  printf("[*] Waiting for client connection...\n");
  
  struct sockaddr_in client_addr;
  socklen_t client_addr_len = sizeof(client_addr);
  client_sock = accept(listen_sock, (struct sockaddr *)&client_addr, &client_addr_len);
  if (client_sock < 0) {
    perror("accept() failed");
    close(listen_sock);
    #ifdef _WIN32
      WSACleanup();
    #endif
    return 1;
  }
  
  printf("[+] Client connected from %s:%d\n\n", inet_ntoa(client_addr.sin_addr),
         ntohs(client_addr.sin_port));
  
  /* ============ STAGE 1: Receive Public Key from Client ============ */
  printf("[STAGE 1] Receiving public key from client...\n");
  
  uint64_t recv_pk_start = get_time_ms();
  if (receive_public_key(client_sock, pk) < 0) {
    fprintf(stderr, "Failed to receive public key\n");
    goto cleanup;
  }
  uint64_t recv_pk_end = get_time_ms();
  
  printf("- Client's Public Key Size: %d bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("[+] Receive PK time: %llu ms\n\n", 
         (unsigned long long)(recv_pk_end - recv_pk_start));
  
  /* ============ STAGE 2: Load Challenge & Send ============ */
  printf("[STAGE 2] Loading challenge and sending to client...\n");
  
  FILE *fin = fopen("test/input.txt", "rb");
  if (!fin) {
    printf("[WARNING] Cannot open test/input.txt, using default challenge\n");
    strcpy((char *)challenge, "This is a test challenge message");
    challenge_len = strlen((char *)challenge);
  } else {
    challenge_len = fread(challenge, 1, BUFFER_SIZE, fin);
    fclose(fin);
  }
  
  printf("- Challenge loaded: %zu bytes\n", challenge_len);
  
  uint64_t send_challenge_start = get_time_ms();
  if (send_challenge(client_sock, challenge, challenge_len) < 0) {
    fprintf(stderr, "Failed to send challenge\n");
    goto cleanup;
  }
  uint64_t send_challenge_end = get_time_ms();
  
  printf("[+] Send challenge time: %llu ms\n\n", 
         (unsigned long long)(send_challenge_end - send_challenge_start));
  
  /* ============ STAGE 3: Receive Signature ============ */
  printf("[STAGE 3] Receiving signature from client...\n");
  
  uint64_t recv_sig_start = get_time_ms();
  if (receive_signature(client_sock, signature, &sig_len) < 0) {
    fprintf(stderr, "Failed to receive signature\n");
    goto cleanup;
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
  
  printf("[+] %s Signature verification successful!\n", 
         verify_result == 0 ? "✓" : "✗");
  
cleanup:
  if (client_sock >= 0) close(client_sock);
  if (listen_sock >= 0) close(listen_sock);
  
  #ifdef _WIN32
    WSACleanup();
  #endif
  
  return verify_result == 0 ? 0 : 1;
}


