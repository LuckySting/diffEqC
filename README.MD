Если что-то поменял в main.c, то
```shell
docker build . -t gccapp
```
Чтобы запустить сбилженный
```shell
docker run -it --rm gccapp
```