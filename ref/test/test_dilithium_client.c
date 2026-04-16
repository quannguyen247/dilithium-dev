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
#define SERVER_IP "127.0.0.1"
#define SERVER_PORT 5000
#define BUFFER_SIZE 8192
#define CHALLENGE_SIZE 1024  /* Size of challenge message */

/* Forward declarations */
uint64_t get_time_ms(void);
int connect_to_server(const char *server_ip, int server_port);
int send_public_key(int sock, const uint8_t *pk);
int receive_challenge(int sock, uint8_t *challenge, size_t *challenge_len, uint64_t *recv_time_ms);
int sign_challenge(const uint8_t *challenge, size_t challenge_len,
                   const uint8_t *sk,
                   uint8_t *signature, size_t *sig_len,
                   uint64_t *sign_time_ms);
int send_signature(int sock, const uint8_t *signature, size_t sig_len);

/* Get current time in milliseconds */
uint64_t get_time_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
}

/* Create TCP socket and connect to server */
int connect_to_server(const char *server_ip, int server_port) {
  int sock;
  struct sockaddr_in server_addr;
  
  printf("[*] Connecting to server at %s:%d\n", server_ip, server_port);
  
  /* Create socket */
  sock = socket(AF_INET, SOCK_STREAM, 0);
  if (sock < 0) {
    perror("socket() failed");
    return -1;
  }
  
  /* Set up server address structure */
  memset(&server_addr, 0, sizeof(server_addr));
  server_addr.sin_family = AF_INET;
  server_addr.sin_port = htons(server_port);
  
  /* Convert IP address */
  if (inet_pton(AF_INET, server_ip, &server_addr.sin_addr) <= 0) {
    perror("inet_pton() failed");
    close(sock);
    return -1;
  }
  
  /* Connect to server */
  if (connect(sock, (struct sockaddr *)&server_addr, sizeof(server_addr)) < 0) {
    perror("connect() failed");
    close(sock);
    return -1;
  }
  
  printf("[+] Successfully connected to server\n");
  return sock;
}

/* Send public key to server */
int send_public_key(int sock, const uint8_t *pk) {
  printf("[*] Sending public key to server (size: %d bytes)...\n", CRYPTO_PUBLICKEYBYTES);
  
  /* Send public key size first */
  uint32_t pk_size = CRYPTO_PUBLICKEYBYTES;
  if (send(sock, (char *)&pk_size, sizeof(pk_size), 0) < 0) {
    perror("send() public key size failed");
    return -1;
  }
  
  /* Send public key */
  ssize_t sent = send(sock, (char *)pk, CRYPTO_PUBLICKEYBYTES, 0);
  if (sent < 0) {
    perror("send() public key failed");
    return -1;
  }
  
  if (sent != CRYPTO_PUBLICKEYBYTES) {
    fprintf(stderr, "Warning: Sent %d bytes, expected %d bytes\n", (int)sent, CRYPTO_PUBLICKEYBYTES);
  }
  
  printf("[+] Public key sent successfully\n");
  return 0;
}

/* Receive challenge from server */
int receive_challenge(int sock, uint8_t *challenge, size_t *challenge_len, uint64_t *recv_time_ms) {
  uint32_t challenge_size;
  
  printf("[*] Waiting for challenge from server...\n");
  
  uint64_t start_time = get_time_ms();
  
  /* Receive challenge size first */
  ssize_t recvd = recv(sock, (char *)&challenge_size, sizeof(challenge_size), 0);
  if (recvd < 0) {
    perror("recv() challenge size failed");
    return -1;
  }
  
  if (recvd == 0) {
    fprintf(stderr, "Connection closed by server\n");
    return -1;
  }
  
  if (challenge_size > BUFFER_SIZE) {
    fprintf(stderr, "Challenge size too large: %u\n", challenge_size);
    return -1;
  }
  
  /* Receive challenge */
  recvd = recv(sock, (char *)challenge, challenge_size, 0);
  if (recvd < 0) {
    perror("recv() challenge failed");
    return -1;
  }
  
  uint64_t end_time = get_time_ms();
  *recv_time_ms = end_time - start_time;
  
  *challenge_len = recvd;
  printf("[+] Challenge received successfully (size: %zu bytes)\n", *challenge_len);
  printf("[+] Challenge reception time: %llu ms\n", (unsigned long long)*recv_time_ms);
  
  return 0;
}

/* Sign the challenge */
int sign_challenge(const uint8_t *challenge, size_t challenge_len,
                   const uint8_t *sk,
                   uint8_t *signature, size_t *sig_len,
                   uint64_t *sign_time_ms) {
  printf("[*] Signing challenge (challenge size: %zu bytes)...\n", challenge_len);
  
  uint64_t start_time = get_time_ms();
  
  /* Sign the challenge */
  int ret = crypto_sign_signature(signature, sig_len,
                                  challenge, challenge_len,
                                  NULL, 0, sk);
  
  uint64_t end_time = get_time_ms();
  *sign_time_ms = end_time - start_time;
  
  if (ret != 0) {
    fprintf(stderr, "crypto_sign_signature() failed with code %d\n", ret);
    return -1;
  }
  
  printf("[+] Challenge signed successfully (signature size: %zu bytes)\n", *sig_len);
  printf("[+] Signing time: %llu ms\n", (unsigned long long)*sign_time_ms);
  
  return 0;
}

/* Send signature to server */
int send_signature(int sock, const uint8_t *signature, size_t sig_len) {
  printf("[*] Sending signature to server (size: %zu bytes)...\n", sig_len);
  
  /* Send signature size first */
  uint32_t size = (uint32_t)sig_len;
  if (send(sock, (char *)&size, sizeof(size), 0) < 0) {
    perror("send() signature size failed");
    return -1;
  }
  
  /* Send signature */
  ssize_t sent = send(sock, (char *)signature, sig_len, 0);
  if (sent < 0) {
    perror("send() signature failed");
    return -1;
  }
  
  if (sent != (ssize_t)sig_len) {
    fprintf(stderr, "Warning: Sent %d bytes, expected %zu bytes\n", (int)sent, sig_len);
  }
  
  printf("[+] Signature sent successfully\n");
  return 0;
}

int main(int argc, char *argv[]) {
  int sock = -1;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t challenge[BUFFER_SIZE];
  size_t challenge_len;
  uint8_t signature[BUFFER_SIZE];
  size_t sig_len;
  
  uint64_t keygen_time_ms;
  uint64_t recv_time_ms;
  uint64_t sign_time_ms;
  
  const char *server_ip = SERVER_IP;
  int server_port = SERVER_PORT;
  
  /* Windows socket initialization */
  #ifdef _WIN32
    WSADATA wsa_data;
    if (WSAStartup(MAKEWORD(2, 2), &wsa_data) != 0) {
      fprintf(stderr, "WSAStartup failed\n");
      return 1;
    }
  #endif
  
  /* Parse command line arguments */
  if (argc > 1) {
    server_ip = argv[1];
  }
  if (argc > 2) {
    server_port = atoi(argv[2]);
  }
  
  printf("\n========== Dilithium Client ==========\n");
  printf("Server: %s:%d\n", server_ip, server_port);
  printf("======================================\n\n");
  
  /* ============ STAGE 1: Generate KeyPair ============ */
  printf("[STAGE 1] Generating keypair...\n");
  printf("- Public Key Size: %d bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("- Secret Key Size: %d bytes\n", CRYPTO_SECRETKEYBYTES);
  
  uint64_t keygen_start = get_time_ms();
  crypto_sign_keypair(pk, sk);
  uint64_t keygen_end = get_time_ms();
  keygen_time_ms = keygen_end - keygen_start;
  
  printf("[+] Keypair generation time: %llu ms\n\n", (unsigned long long)keygen_time_ms);
  
  /* ============ STAGE 2: Connect to Server & Send Public Key ============ */
  printf("[STAGE 2] Connecting to server and sending public key...\n");
  
  sock = connect_to_server(server_ip, server_port);
  if (sock < 0) {
    fprintf(stderr, "Failed to connect to server\n");
    #ifdef _WIN32
      WSACleanup();
    #endif
    exit(EXIT_FAILURE);
  }
  
  if (send_public_key(sock, pk) < 0) {
    fprintf(stderr, "Failed to send public key\n");
    close(sock);
    #ifdef _WIN32
      WSACleanup();
    #endif
    exit(EXIT_FAILURE);
  }
  
  printf("\n");
  
  /* ============ STAGE 3: Wait for Challenge ============ */
  printf("[STAGE 3] Waiting for challenge from server...\n");
  
  if (receive_challenge(sock, challenge, &challenge_len, &recv_time_ms) < 0) {
    fprintf(stderr, "Failed to receive challenge\n");
    close(sock);
    #ifdef _WIN32
      WSACleanup();
    #endif
    exit(EXIT_FAILURE);
  }
  
  printf("\n");
  
  /* ============ STAGE 4: Sign the Challenge ============ */
  printf("[STAGE 4] Signing the challenge...\n");
  printf("- Signature Size: %d bytes\n", CRYPTO_BYTES);
  
  if (sign_challenge(challenge, challenge_len, sk, signature, &sig_len, &sign_time_ms) < 0) {
    fprintf(stderr, "Failed to sign challenge\n");
    close(sock);
    #ifdef _WIN32
      WSACleanup();
    #endif
    exit(EXIT_FAILURE);
  }
  
  printf("\n");
  
  /* ============ STAGE 5: Send Signature to Server ============ */
  printf("[STAGE 5] Sending signature to server...\n");
  
  if (send_signature(sock, signature, sig_len) < 0) {
    fprintf(stderr, "Failed to send signature\n");
    close(sock);
    #ifdef _WIN32
      WSACleanup();
    #endif
    exit(EXIT_FAILURE);
  }
  
  close(sock);
  
  printf("\n");
  
  /* ============ RESULTS ============ */
  printf("===================================\n");
  printf("[TIMING SUMMARY]\n");
  printf("===================================\n");
  printf("KeyGen Time:               %llu ms\n", (unsigned long long)keygen_time_ms);
  printf("Challenge Reception Time:  %llu ms\n", (unsigned long long)recv_time_ms);
  printf("Signing Time:              %llu ms\n", (unsigned long long)sign_time_ms);
  printf("Total Time (Gen + Sign):   %llu ms\n", 
         (unsigned long long)(keygen_time_ms + sign_time_ms));
  printf("===================================\n\n");
  
  printf("[KEY INFORMATION]\n");
  printf("- Public Key Size:  %d bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("- Secret Key Size:  %d bytes\n", CRYPTO_SECRETKEYBYTES);
  printf("- Signature Size:   %zu bytes\n", sig_len);
  printf("- Challenge Size:   %zu bytes\n", challenge_len);
  printf("===================================\n\n");
  
  printf("[+] Client completed successfully!\n");
  
  /* Windows socket cleanup */
  #ifdef _WIN32
    WSACleanup();
  #endif

  return 0;
}
