#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

#include "../randombytes.h"
#include "../sign.h"

#define SERVER_PORT 5000
#define BUFFER_SIZE 10000
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

int send_packet(int sock, uint8_t type, const uint8_t *data, uint32_t data_len) {
    PacketHeader header;
    header.magic = htons(MAGIC_NUMBER);
    header.type = type;
    header.length = htonl(data_len);

    if (send(sock, (char *)&header, sizeof(header), 0) < 0) return -1;
    if (send(sock, (char *)data, data_len, 0) < 0) return -1;
    return 0;
}

int receive_packet(int sock, uint8_t *buffer, uint32_t *len) {
    PacketHeader header;
    if (recv(sock, (char *)&header, sizeof(header), 0) <= 0) return -1;

    uint32_t payload_len = ntohl(header.length);
    uint32_t total_received = 0;
    while (total_received < payload_len) {
        int n = recv(sock, (char *)buffer + total_received, payload_len - total_received, 0);
        if (n <= 0) break;
        total_received += n;
    }
    *len = total_received;
    return header.type;
}

int main(int argc, char *argv[]) {
    int sock;
    struct sockaddr_in server_addr;
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t challenge[BUFFER_SIZE];
    uint32_t challenge_len;
    uint8_t signature[CRYPTO_BYTES];
    size_t sig_len;

    char *ip = (argc > 1) ? argv[1] : "127.0.0.1";

    printf("[*] Generating Dilithium keys...\n");
    crypto_sign_keypair(pk, sk);

    sock = socket(AF_INET, SOCK_STREAM, 0);
    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(SERVER_PORT);
    inet_pton(AF_INET, ip, &server_addr.sin_addr);

    if (connect(sock, (struct sockaddr *)&server_addr, sizeof(server_addr)) < 0) {
        perror("Connect failed");
        return 1;
    }
    printf("[+] Connected to %s\n", ip);

    printf("[*] Sending Public Key...\n");
    send_packet(sock, TYPE_PK, pk, CRYPTO_PUBLICKEYBYTES);

    printf("[*] Waiting for Challenge from Python Server...\n");
    receive_packet(sock, challenge, &challenge_len);
    printf("[+] Received challenge: %u bytes\n", challenge_len);

    printf("[*] Signing challenge...\n");
    crypto_sign_signature(signature, &sig_len, challenge, challenge_len, NULL, 0, sk);

    printf("[*] Sending Signature...\n");
    send_packet(sock, TYPE_SIG, signature, (uint32_t)sig_len);

    printf("[DONE] Dilithium process finished successfully.\n");
    close(sock);
    return 0;
}