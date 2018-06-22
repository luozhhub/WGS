#/bin/bash
# 以 shell 方式启动容器，环境如同 Linux 命令行

# 将 vnpy 指定到
docker run --name mysql  -p 5900:5900 -v `pwd`:/root/UMLS/UMLS-2016AB/2016AB -v /simm/home/luozhihui/try/umls:/var/lib/mysql/umls_mysql/ -it luozh/umls_mysql:latest
