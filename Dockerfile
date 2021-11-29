FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y libblas-dev liblapacke-dev liblapack-dev gfortran
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp
RUN gcc main.c -lm -llapack -llapacke -lblas -lgfortran -o my_prog
CMD ["./my_prog"]