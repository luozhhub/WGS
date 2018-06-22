FROM ubuntu:16.04

# 时区设置
ENV TZ=Asia/Shanghai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# 中文支持
RUN apt-get clean && apt-get update
RUN apt-get install locales
RUN locale-gen zh_CN.UTF-8 \
    && DEBIAN_FRONTEND=noninteractive dpkg-reconfigure locales
RUN locale-gen zh_CN.UTF-8
ENV LANG zh_CN.UTF-8
ENV LANGUAGE zh_CN:zh
ENV LC_ALL zh_CN.UTF-8


## 替换 apt 源为阿里云，在本地构建镜像时，取消注释，使用阿里云的 apt 源
#RUN echo "" > /etc/apt/sources.list \
#    && echo "deb http://mirrors.aliyun.com/ubuntu/ xenial main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb http://mirrors.aliyun.com/ubuntu/ xenial-backports main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb http://mirrors.aliyun.com/ubuntu/ xenial-proposed main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb http://mirrors.aliyun.com/ubuntu/ xenial-security main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb http://mirrors.aliyun.com/ubuntu/ xenial-updates main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb-src http://mirrors.aliyun.com/ubuntu/ xenial main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb-src http://mirrors.aliyun.com/ubuntu/ xenial-backports main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb-src http://mirrors.aliyun.com/ubuntu/ xenial-proposed main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb-src http://mirrors.aliyun.com/ubuntu/ xenial-security main multiverse restricted universe" >> /etc/apt/sources.list \
#    && echo "deb-src http://mirrors.aliyun.com/ubuntu/ xenial-updates main multiverse restricted universe" >> /etc/apt/sources.list

# 更新源列表
RUN apt-get clean \
    && apt-get update

# 安装编译环境
