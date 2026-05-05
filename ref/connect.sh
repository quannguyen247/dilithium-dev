#!/bin/bash

TARGET_IP="192.168.4.85"
CLIENT_BIN="./test/test_dilithium_client2"
CONCURRENT_SESSIONS=10

echo "--- START TESTING ON $TARGET_IP ---"

while true; do
    echo "Đang khởi tạo $CONCURRENT_SESSIONS kết nối đồng thời..."

    for i in $(seq 1 $CONCURRENT_SESSIONS); do
        # Chạy client dưới nền và ẩn output
        $CLIENT_BIN "$TARGET_IP" > /dev/null 2>&1 &
    done

    # Chờ tất cả các tiến trình con trong đợt này hoàn thành
    wait

    echo "=> Đã xử lý xong một đợt!"
    echo "-----------------------------------------------"
done