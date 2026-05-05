#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <limits.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

#include "../randombytes.h"
#include "../sign.h"

#define SERVER_PORT 5000
#define BUFFER_SIZE 8192

static int send_all(int sock, const uint8_t *buf, size_t len) {
    size_t total = 0;
    while (total < len) {
        ssize_t sent = send(sock, (const char *)buf + total, (int)(len - total), 0);
        if (sent < 0) {
            if (errno == EINTR) {
                continue;
            }
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
            if (errno == EINTR) {
                continue;
            }
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

static int recv_blob(int sock, uint8_t *buffer, uint32_t buffer_size, uint32_t *out_len) {
    uint32_t len_net = 0;
    if (recv_all(sock, (uint8_t *)&len_net, sizeof(len_net)) < 0) {
        return -1;
    }

    uint32_t payload_len = ntohl(len_net);
    if (payload_len > buffer_size) {
        errno = EMSGSIZE;
        return -1;
    }

    if (payload_len > 0 && recv_all(sock, buffer, payload_len) < 0) {
        return -1;
    }

    *out_len = payload_len;
    return 0;
}

int main(int argc, char *argv[]) {
    int sock = -1;
    struct sockaddr_in server_addr;
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t challenge[BUFFER_SIZE];
    uint32_t challenge_len = 0;
    uint8_t signature[CRYPTO_BYTES];
    size_t sig_len = 0;

    const char *ip = (argc > 1) ? argv[1] : "127.0.0.1";

    printf("[*] Generating Dilithium keys...\n");
    if (crypto_sign_keypair(pk, sk) != 0) {
        fprintf(stderr, "Key generation failed\n");
        return 1;
    }

    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
        perror("socket failed");
        return 1;
    }

    memset(&server_addr, 0, sizeof(server_addr));
    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(SERVER_PORT);
    if (inet_pton(AF_INET, ip, &server_addr.sin_addr) != 1) {
        fprintf(stderr, "Invalid server IP: %s\n", ip);
        close(sock);
        return 1;
    }

    if (connect(sock, (struct sockaddr *)&server_addr, sizeof(server_addr)) < 0) {
        perror("connect failed");
        close(sock);
        return 1;
    }
    printf("[+] Connected to %s\n", ip);

    printf("[*] Sending public key...\n");
    if (send_blob(sock, pk, CRYPTO_PUBLICKEYBYTES) < 0) {
        perror("send() public key failed");
        close(sock);
        return 1;
    }

    printf("[*] Waiting for challenge from server...\n");
    if (recv_blob(sock, challenge, BUFFER_SIZE, &challenge_len) < 0) {
        perror("recv() challenge failed");
        close(sock);
        return 1;
    }
    printf("[+] Received challenge: %u bytes\n", challenge_len);

    printf("[*] Signing challenge...\n");
    if (crypto_sign_signature(signature, &sig_len, challenge, (size_t)challenge_len, NULL, 0, sk) != 0) {
        fprintf(stderr, "Signature failed\n");
        close(sock);
        return 1;
    }

    if (sig_len > UINT32_MAX) {
        fprintf(stderr, "Signature too large: %zu bytes\n", sig_len);
        close(sock);
        return 1;
    }

    printf("[*] Sending signature...\n");
    if (send_blob(sock, signature, (uint32_t)sig_len) < 0) {
        perror("send() signature failed");
        close(sock);
        return 1;
    }

    printf("[DONE] Dilithium process finished successfully.\n");
    close(sock);
    return 0;
}