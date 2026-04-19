#!/bin/sh
IP="192.168.4.85"

echo "BẮT ĐẦU BAN TUM LUM TỪ CLIENT..."

while true
do
    # Bắn 10 client cùng một lúc vào background (dùng dấu &)
    for i in 1 2 3 4 5 6 7 8 9 10
    do
        ./test/test_dilithium_client2 $IP > /dev/null 2>&1 &
    done
    
    # Đợi 10 client này chạy xong thì mới bắn tiếp đợt mới
    wait
    echo "Đã xử lý xong 1 đợt 10 kết nối đồng thời!"
done
