#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h> /* Thư viện đa luồng */

/* Platform-specific socket headers - must come before ../sign.h to avoid m>
#ifdef _WIN32
  #include <winsock2.h>
  #include <ws2tcpip.h>
  #define close(sock) closesocket(sock)
  #define ssize_t int
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

/* Dilithium headers */
#include "../randombytes.h"
#include "../sign.h"

/* Configuration */
#define SERVER_PORT 5000
#define BUFFER_SIZE 8192
#define MLEN 1200  /* Challenge size */

#define MAGIC_NUMBER 0xABCD
#define TYPE_PK    1
#define TYPE_CHAL  2
#define TYPE_SIG   3

#pragma pack(push, 1)
typedef struct {
    uint16_t magic;
    uint8_t  type;
    uint32_t length;
} PacketHeader;
#pragma pack(pop)

/* Cấu trúc dữ liệu truyền vào Thread */
typedef struct {
    int socket;
    struct sockaddr_in address;
} client_info_t;

/* Forward declarations */
uint64_t get_time_ms(void);
int receive_public_key(int sock, uint8_t *pk);
int send_challenge(int sock, const uint8_t *challenge, size_t challenge_len>
int receive_signature(int sock, uint8_t *signature, size_t *sig_len);
/* Get current time in milliseconds */
uint64_t get_time_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
}

/* Receive public key from client */
int receive_public_key(int sock, uint8_t *pk) {
  PacketHeader header;
  if (recv(sock, (char *)&header, sizeof(header), 0) <= 0) return -1;

  if (ntohs(header.magic) != MAGIC_NUMBER) {
      fprintf(stderr, "Invalid Magic Number!\n");
      return -1;
  }

  uint32_t pk_size = ntohl(header.length);
  if (pk_size != (uint32_t)CRYPTO_PUBLICKEYBYTES) {
    fprintf(stderr, "Invalid public key size: %u (expected %d)\n", pk_size,>
    return -1;
  }

  uint32_t total = 0;
  while (total < pk_size) {
      ssize_t n = recv(sock, (char *)pk + total, pk_size - total, 0);
      if (n <= 0) return -1;
      total += n;
  }
  return 0;
}

/* Send challenge message to client */
int send_challenge(int sock, const uint8_t *challenge, size_t challenge_len>
  PacketHeader header;
header.magic = htons(MAGIC_NUMBER);
  header.type = TYPE_CHAL;
  header.length = htonl((uint32_t)challenge_len);

  if (send(sock, (char *)&header, sizeof(header), 0) < 0) return -1;
  if (send(sock, (char *)challenge, challenge_len, 0) < 0) return -1;
  return 0;
}

/* Receive signature from client */
int receive_signature(int sock, uint8_t *signature, size_t *sig_len) {
  PacketHeader header;
  if (recv(sock, (char *)&header, sizeof(header), 0) <= 0) return -1;

  uint32_t size = ntohl(header.length);
  if (size > BUFFER_SIZE) {
    fprintf(stderr, "Signature size too large: %u\n", size);
    return -1;
  }

  uint32_t total = 0;
  while (total < size) {
      ssize_t n = recv(sock, (char *)signature + total, size - total, 0);
      if (n <= 0) return -1;
      total += n;
  }

  *sig_len = (size_t)total;
  return 0;
}

/* ========================================================= */
/* LUỒNG XỬ LÝ RIÊNG CHO TỪNG BOARD (WORKER THREAD)          */
/* ========================================================= */
void *handle_client(void *arg) {
    client_info_t *info = (client_info_t *)arg;
int client_sock = info->socket;
    char client_ip[INET_ADDRSTRLEN];
    inet_ntop(AF_INET, &info->address.sin_addr, client_ip, INET_ADDRSTRLEN);

    printf("\n[+] Client connected from %s:%d\n\n", client_ip, ntohs(info->>

    uint64_t total_start = get_time_ms();
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t challenge[BUFFER_SIZE];
    size_t challenge_len;
    uint8_t signature[BUFFER_SIZE];
    size_t sig_len;

    /* ============ STAGE 1 ============ */
    printf("[%s - STAGE 1] Receiving public key...\n", client_ip);
    uint64_t recv_pk_start = get_time_ms();
    if (receive_public_key(client_sock, pk) < 0) goto cleanup;
    uint64_t recv_pk_end = get_time_ms();
    printf("[%s] Receive PK time: %llu ms\n\n", client_ip, (unsigned long l>

    /* --- ĐÓNG BĂNG 15 GIÂY ĐỂ BẠN KỊP BẬT CÁC BOARD KHÁC --- */
    printf(">>>>> [%s] ĐANG CHỜ 15 GIÂY... HÃY BẬT CÁC BOARD KHÁC LÊN! <<<<>

    /* ------------------------------------------------------- */

    /* ============ STAGE 2 ============ */
    printf("[%s - STAGE 2] Loading and sending challenge...\n", client_ip);
    FILE *fin = fopen("test/input.txt", "rb");
    if (!fin) {
        strcpy((char *)challenge, "This is a test challenge message");
        challenge_len = strlen((char *)challenge);
    } else {
        challenge_len = fread(challenge, 1, BUFFER_SIZE, fin);
        fclose(fin);
    }
 uint64_t send_challenge_start = get_time_ms();
    if (send_challenge(client_sock, challenge, challenge_len) < 0) goto cle>
    uint64_t send_challenge_end = get_time_ms();
    printf("[%s] Send challenge time: %llu ms\n\n", client_ip, (unsigned lo>

    /* ============ STAGE 3 ============ */
    printf("[%s - STAGE 3] Receiving signature...\n", client_ip);
    uint64_t recv_sig_start = get_time_ms();
    if (receive_signature(client_sock, signature, &sig_len) < 0) goto clean>
    uint64_t recv_sig_end = get_time_ms();
    printf("[%s] Receive signature time: %llu ms\n\n", client_ip, (unsigned>

    /* ============ STAGE 4 ============ */
    printf("[%s - STAGE 4] Verifying signature...\n", client_ip);
    uint64_t verify_start = get_time_ms();
    int verify_result = crypto_sign_verify(signature, sig_len, challenge, c>
    uint64_t verify_end = get_time_ms();

    uint64_t total_end = get_time_ms();

    /* ============ TIMING SUMMARY ============ */
    printf("===================================\n");
    printf("[%s - TIMING SUMMARY]\n", client_ip);
    printf("Receive PK Time:           %llu ms\n", (unsigned long long)(rec>
    printf("Send Challenge Time:       %llu ms\n", (unsigned long long)(sen>
    printf("Receive Signature Time:    %llu ms\n", (unsigned long long)(rec>
    printf("Verification Time:         %llu ms\n", (unsigned long long)(ver>
    printf("Total Time (from start):   %llu ms\n", (unsigned long long)(tot>
    printf("===================================\n");
    printf("[+] %s Signature verification successful for %s!\n\n", verify_r>

cleanup:
    close(client_sock);
    free(info); // Giải phóng RAM cho thread
    return NULL;
}


/* ========================================================= */
/* LUỒNG CHÍNH (MAIN THREAD) - CHỈ LÀM NHIỆM VỤ ĐÓN KHÁCH    */
/* ========================================================= */
int main(void) {
  int listen_sock = -1;

  printf("\n========== Multi-threaded Dilithium Server ==========\n");
  printf("Listening on port %d\n", SERVER_PORT);
  printf("=====================================================\n\n");

  listen_sock = socket(AF_INET, SOCK_STREAM, 0);
  int reuse = 1;
  setsockopt(listen_sock, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(>

  struct sockaddr_in server_addr;
  memset(&server_addr, 0, sizeof(server_addr));
  server_addr.sin_family = AF_INET;
  server_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  server_addr.sin_port = htons(SERVER_PORT);

  bind(listen_sock, (struct sockaddr *)&server_addr, sizeof(server_addr));
  listen(listen_sock, 10); // Cho phép xếp hàng 10 kết nối

  printf("[*] Server is running and waiting for clients...\n");

  /* VÒNG LẶP VÔ TẬN - KHÔNG BAO GIỜ TẮT SERVER */
  while (1) {
      struct sockaddr_in client_addr;
      socklen_t client_addr_len = sizeof(client_addr);
/ Luồng chính đứng chờ ở đây
      int client_sock = accept(listen_sock, (struct sockaddr *)&client_addr>
      if (client_sock < 0) continue;

      // Chuẩn bị dữ liệu để quăng cho Luồng phụ (Worker)
      client_info_t *info = malloc(sizeof(client_info_t));
      info->socket = client_sock;
      info->address = client_addr;

      // Tạo Luồng phụ để xử lý Board này
      pthread_t thread_id;
      if (pthread_create(&thread_id, NULL, handle_client, (void *)info) != >
          perror("Failed to create thread");
          close(client_sock);
          free(info);
      } else {
          // Tách luồng để nó tự lo liệu vòng đời của mình
          pthread_detach(thread_id);
      }
  }

  close(listen_sock);
  return 0;
}