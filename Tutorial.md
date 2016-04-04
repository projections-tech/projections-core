#Projections Tutorial

This tutorial covers Projections meta filesystem basic concepts and use.

##Installation
###Requirements
This demo was tested on Ubuntu 15.04 x64 (in order to run Docker), newer versions are compatible too.
###Installation instructions
1) In this step Docker is set up, in order to run projections database container, if Docker is already installed on 
system user may skip this step. In order to perform Docker installation, firstly user should update APT package index:

```
$ sudo apt-get update
```

Then install Docker itself and start installed service:

```
$ sudo apt-get install docker-engine && service docker start
```

And add current user (substitute {user_name} with current user name) to docker usergroup in order to run it without sudo:

```
$ sudo usermod -aG docker {user_name}
```

After this command finishes it is required to perfrom log out and log back in, in order to ensure that user 
is running with the correct permissions.

2) Clone projections demo using git clone. If git is not installed on host, user need to run following command:

```
$ sudo apt-get install git
```

And then perform git setup as follows, substituting name and email with current user name:

```
$ git config --global user.name {user_name}
```

And email:

```
$ git config --global user.email {user_email}
```

After this step user need to clone demo using git "clone" command:

```
$ git clone https://github.com/projections-tech/projections-core.git
```

If user works in, for example "home" directory, when demo will be cloned into "~/projections-core" directory. 
In order to perform next operations user need`s to change directory into cloned directory:

```
$ cd ~/projections-core
```

3) In this step projection database image created and it`s container started. To build projections database image 
following command is used:

```
$ docker build -t projections_database -f database.dockerfile .
```

In order to start projections database container following command is used:

```
$ docker run -d -p 48620:5432 --name projections_database projections_database
```

Here, parameter "-p 48621:5432" binds database container port 5432 to port 48261 on host. 

*WARNING! It is assumed that port 48261 is not busy on host, otherwise user must specify desired port manually.*

4) Final step is to create python virtual environment in which demo will run and start it. If user does`t have 
virtualenv installed, installation performed using following command:

```
$ sudo apt-get install virtualenv
```

When virtualenv installation is finished, demo environment is started using following command:

```
$ virtualenv -p python3 --no-site-packages --distribute projections_demo  && source projections_demo/bin/activate && pip install -r requirements.txt
```

When this command is finished the demo setup is complete. In order to exit demo virtual environment, user need to type:

```
$ deactivate
```

In order to reenter environment following command is entered from projections demo directory:

```
$ source projections_demo/bin/activate
```

Database container is stopped using Docker command:

```
$ docker stop projections_database
```

And removed using following command:

```
$ docker rm projections_database
```

In order to remove projections database image from system user needs to type:

```
$ docker rmi projections_database
```

##Basic projections usage
###Starting projections daemon
***Note!*** In following examples it is assumed that projections demo where installed into ~/projections-core.

On it`s most basic level, projection - is filesystem representation of various heterogeneous resources.

We begin our tutorial by starting projections daemon. In projections meta filesystem, daemon is entity that manages 
relationships between projectors - individual processes, which create resource projections and projection database - 
projections persistance storage.

Work with the daemon is quite simple, we can start it using following command:

```
$ ./projections_daemon.py -start
```

