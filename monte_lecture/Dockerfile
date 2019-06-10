FROM ubuntu

WORKDIR /home/files
COPY ./files /home/files
RUN apt update -y && apt upgrade -y && apt install -y gnuplot

CMD ["gnuplot", "plot-cos.plt"]
