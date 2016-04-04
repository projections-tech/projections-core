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

##Projections usage
###Starting projections daemon
***Note!*** *In following examples it is assumed that projections demo where installed into ~/projections-core.*

On it`s most basic level, projection - is filesystem representation of various heterogeneous resources.

We begin our tutorial by starting projections daemon. In projections meta filesystem, daemon is entity that manages 
relationships between projectors - individual processes, which create resource projections and projection database - 
projections persistent storage.

Work with the daemon is quite simple, we can start it using following command:

```
$ ./projections_daemon.py -start
```

Daemon will prompt:

```
Starting projections daemon!
```

Similarly we can stop daemon when we finish work with projections using this command:

```
$ ./projections_daemon.py -stop
```

Daemon is intended to work indefinitely, managing users projections and it`s worth to mention that only one daemon process is allowed per host.

Now, with projections daemon started, let`s jump right to action.

###Basic projections usage
As stated earlier projection is representation of various resources. These resources may include, for example, contents 
of some remote filesystem, data from scientific resources like NCBI Genbank, results of sequencing run and even local 
filesystem contents.
####Creating projection
Let`s for example, project genbank query as set of filesystem objects. To do this, we issue following command to 
projections daemon command line interface:

```
$ ./projections_cli.py project --projection_name=genbank_example --mount_point=mount --driver=genbank_driver --prototype=examples/genbank_example_1_config.yaml
```

Let`s describe command arguments:

- First argument, "project" is self evident, we send command to daemon to create projection.
- "--projection_name" - second argument states the name of created projection, using which we can differentiate our 
projection from others. It is important to mention that projection name **must** be unique. Shorthand for this 
command is "-n".
- "--mount_point" - this argument specifies mount point on which projections will be mounted. We will see created projections 
there. Mount point **must** also be unique for each created projection. Shorthand for this command is "-m".
- "--driver" - this argument specifies driver, which will be used to interact with resource that is projected. 
We choose to project Genbank so choice is obvious. Shorthand for this command is "-d".
- "--prototype" - this argument specifies projection prototype. In order to encompass such a variety of resources, 
projections use small text files, which contain resource logical structure - the *prototype* of resource to be projected. 
We will return to this topic shortly, for now, let`s just say that we want to project 5 first matches for query to 
Genbank searching for "escherichia". Shorthand for this command is "-p". 

When projection is created, daemon will prompt:

```
Projection "genbank_example" created and started!
```

Now, using command like "ls -R mount/", or simply moving to directory using nautilus desktop manager we will see following result:

```
mount/:
escherichia[orgn]_1013076717  escherichia[orgn]_1013076742  escherichia[orgn]_1013076809  escherichia[orgn]_1013077104  escherichia[orgn]_1013077116
```

```
mount/escherichia[orgn]_1013076717:
1013076717.fasta  1013076717.gb
```

```
mount/escherichia[orgn]_1013076742:
1013076742.fasta  1013076742.gb
```

```
mount/escherichia[orgn]_1013076809:
1013076809.fasta  1013076809.gb
```

```
mount/escherichia[orgn]_1013077104:
1013077104.fasta  1013077104.gb
```

```
mount/escherichia[orgn]_1013077116:
1013077116.fasta  1013077116.gb
```

Query projection is now created, now we can open files in projections subdirectories using our text editor of choice or perform data analysis.
####Listing projections
Following command is used to list available projections:

```
./projections_cli.py ps
```

And will return list of avaliable projections: 

```
Projection id: 1 Projection name: genbank_example Mount point: mount Driver: genbank_driver Projector PID: 22127
```

Here we can see that this is first projection created in database, hence id is 1, projection name is set as we entered 
earlier, projection mount point is mount/ and driver is genbank driver. Projector PID is id of projector process which 
currently performs projection in directory mount, it will differ from example and set to None when we stop projection.

####Stopping and starting projections
After we checked that projections is fully created and works, let`s perform projection stop using this command:

```
$ ./projection_cli.py stop -n genbank_example 
```

This will result in daemon response:

```
Stopped projection: "genbank_example".
```

Let`s check that projection is already stopped using "ps" command:

```
Projection id: 1 Projection name: genbank_example Mount point: mount Driver: genbank_driver Projector PID: None
```

As we can see, Projector PID is set to None, since no Projector process is running. Projection itself is stored in projection database. 

If we want to restart projection, we need to type:
 
```
$ ./projection_cli.py start -n genbank_example 
```

Which will restart projection.

####Removing projections
As we move to next tutorial part, projection "genbank_example" is no longer needed, we can remove it using "rm" command:

```
./projections_cli.py rm -n genbank_example
```

Which will stop and remove "genbank_example" projection. We also can check this using "ps" function.