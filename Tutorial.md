# Projections Tutorial

This tutorial covers Projections meta filesystem basic concepts and use.

## Installation
### Requirements
This demo was tested on Ubuntu 15.04 x64 (in order to run Docker), newer versions are compatible too.
### Installation instructions
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

If user works in, for example, "home" directory, when demo will be cloned into "~/projections-core" directory. 
In order to perform next operations user have to change directory into cloned directory:

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

Here, parameter "-p 48620:5432" binds database container port 5432 to port 48260 on host. 

*WARNING! It is assumed that port 48260 is not busy on host, otherwise user must specify desired port manually.*

4) Final step is to create python virtual environment in which demo will run and start it. If user doesn`t have 
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

In order to reenter environment following command is entered from projections-core directory:

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

## Projections usage
### Starting projections daemon
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

### Basic projections usage
As stated earlier projection is representation of various resources. These resources may include, for example, contents 
of some remote filesystem, data from scientific resources like NCBI Genbank, results of sequencing run and even local 
filesystem contents.
#### Creating projection
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

mount/escherichia[orgn]_1013076717:
1013076717.fasta  1013076717.gb

mount/escherichia[orgn]_1013076742:
1013076742.fasta  1013076742.gb

mount/escherichia[orgn]_1013076809:
1013076809.fasta  1013076809.gb

mount/escherichia[orgn]_1013077104:
1013077104.fasta  1013077104.gb

mount/escherichia[orgn]_1013077116:
1013077116.fasta  1013077116.gb
```

Query projection is now created, now we can open files in projections subdirectories using our text editor of choice or perform some data analysis.
#### Listing projections
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

#### Stopping and starting projections
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

#### Removing projections
As we move to next tutorial part, projection "genbank_example" is no longer needed, we can remove it using "rm" command:

```
./projections_cli.py rm -n genbank_example
```

Which will stop and remove "genbank_example" projection. We also can check this using "ps" function.

### Prototypes structure

As were briefly mentioned above, projections are created using prototype files, which describe projection logical 
structure. Projection itself is a tree, with nodes corresponding to projected resource elements, like query result was
folder in example above. To create projection nodes tree we define projection prototypes tree - nodes of these tree
hold rules which define how to create current projection node basing on context of upstream projection nodes.

All prototypes share this basic structure:

```
root_projection_uri -> root -> prototype_definition
```

Where: 

- "root_projection_uri" - is root URI for whole projection. In genbank example this was "escherichia" search query. This
uri creates context for all downstream nodes. In some sense this is projection created without prototype by user.
- "root" - root projection creates first level of a projections tree using rules stored in it\`s prototypes. URI\`s of 
these new nodes serve as context for lower level nodes and so on to next prototype definition.

Each prototype have it`s own environment which consist of three fields:

- environment - metadata associated with upper level projection, using code stored in prototype projection URI as resolved from this field. 
- content - metadata of projection at present URI, from this environment prototype resolves projection name. Field "content_uri" is added to this metadata.
- context - metadata of all upstream nodes.

Each prototype in prototypes tree have three fields that describe it: URI, Name and Meta_links:

-"URI" - this field holds code, using which URI of projection will be resolved from prototype`s environment.

-"Name" - this field holds code, using wich Name ofprojection is resolved, according to it`s content.

-"Meta_links" - this field holds dictionary of link names and link codes, which define data-metadata binding in projection. This action is perfromed after projection tree is complete.

Let`s look at actual prototype file, and how it is unfolded into projections step by step 
(file "genbank_example_config.yaml" in "examples/" dir, we read it from bottom to top). Pototype file is written in 
[YAML format](http://www.yaml.org/spec/1.2/spec.html), URI an Name microcodes use 
[ObjectPath query language](http://objectpath.org/) and SQL is used to bind metadata to data.

Before prototypes definition comes header, header ends with root_resource_uri from which projection will bes started, we can try any root uri which acts like NCBI search query (for example search for 'Homo Sapiens[orgn] AND Tumor'):

``` YAML
# Prototype microcode dialects are defined here.
dialect: ObjectPath, SQL
# This is uri of resource which will be projected.
resource_uri: http://eutils.ncbi.nlm.nih.gov
# This is a path to projection drvier configuration file, it may content access credentials, driver-specific settings etc.
driver_config_path: 'drivers/driver_configurations/genbank_config.yaml'

# This is root projection URI, it consists of three parts:
# 1) 'search_query:' suffix, which tells driver to process query as search query
# 2) Search query as accepted by NCBI. All query parameters and operators as described here:
# http://www.ncbi.nlm.nih.gov/books/NBK49540/
# 3) Number of results to return, if no number specified driver will return one query result
root_projection_uri: 'search_query:escherichia[orgn]:5'
```

At first comes root prototype:

``` YAML
# This is root projection node definition. By convention, it must be named "root". By resolving it`s URI code we create
# list of URI`s for associated Genbank entries.
root:
  # Type of root node. Directory contains other projections.
  type: directory
  # This is a dictionary of children prototypes for current prototype.
  children:
    # Here we add an aliases for lower level prototypes.
    gb_file_prototype: *gb_file_prototype
    fasta_file_prototype: *fasta_file_prototype
  # This field encodes how node name is defined from node "content" field we join query name and URI of query.
  name: " join([$.environment.TranslationSet[0].From, $.content_uri], '_') "
  # This field defines how URI is assigned to projections.
  uri: " $.IdList"
```

This prototype unfolds into directories which correspond to search query.

Secondly we define GB file prototype, notice how prototype name contains anchor "&gb_file_prototype", it was used before 
in root prototype:

``` YAML
# This is GB file prototype, it`s link created by resoloving environment which is avaliable by the URI resolved in
# upper level. File extension is appended to complete prototype URI, driver will resolve it accordingly. Name of prototype
# is resolved using prototype context, here it is setted as prototype uri, using field avaliable for every prototype -
# "content_uri"
gb_file_prototype: &gb_file_prototype
  type: file
  # File prototype cannot have children prototypes by definition.
  children: None
  name: " $.content_uri "
  uri: " join([$.IdList[0], '.gb'], '')"
```

Third prototype is FASTA file prototype:

``` YAML
# This is FASTA file prototype, corresponding for result query, it`s created similarly to GB file prototype, by adding
# ".fasta" extension to prototype URI. This prototype differs from GB prototype by extension, added to link.
fasta_file_prototype: &fasta_file_prototype
  type: file
  children: None
  name: " $.content_uri "
  # Here we define prototype metadata link. Metadata is defined using SQL microcode and is assigned to projection after
  # projection tree is completed. Code describes single node on tree which abides query conditions.
  # Here, code states that metadata node is a node where node name is equal to current projection name and it`s extension is
  # ".gb". This effectively binds GB files to their corresponding FASTA files.
  meta_link:
  - " regexp_replace(nodes.node_name, '.gb', '') = regexp_replace(current_node.node_name, '.fasta', '') "
  uri: " join([$.IdList[0], '.fasta'], '')"
```

You may notice, that presented prototype order is reversed, this is due to alias cannot be created before anchor. Root
node have no anchors attached, and is in lowest position in file.

#### Selectors in prototype definitions
During projection creation we are able to filter which resources to project using selectors in microcode. 
To illustrate selectors, let`s start another projection, this time we will project local filesystem contents:

```
$ ./projections_cli.py project -n fs_example_1 -m mount -p examples/fs_example_1_config.yaml -d fs_driver
```

Now we should examine projection configuration file "examples/fs_example_1_config.yaml". Root uri for this projection 
is "tests/test_dir" - this URI corresponds to test directory in demo folder, we can check if projection is created 
correctly comparing result with it. Root projection resembles Genbank root projection, and contains 4 children prototypes, 
one of them is directory prototype "samples", take a closer look of it`s URI:

``` YAML
uri: " $.children[@.type is 'dir'].resource_uri "
```

URI object path resolution code means following:"For node in root`s children ("$.children part"), select node`s resource URI (".resource_uri" part) 
if node`s type is 'dir'(code in square brackets)". This was example of selector to filter root`s node children according their environment 
fields values. Driver fetches environment according to root URI for "samples" directory, fetches content from URI resolved 
from environment using ObjectPath code using this content "samples" name is set. Children nodes for "samples" will have 
their environment uri as "samples" resolved URI.

Root prototype have other prototypes, we will look at fasta_prototype, other files selectors are similiar:

```
uri: " $.children[@.type is 'file' and @.extension is '.fasta'].resource_uri "
```

This code means following: "For node in sample prototype`s children, select URI of node which type is 'file' and node extension 
is .fasta". ObjectPath also allows use of other helper functions in selectors, reader is encouraged to read it`s 
[documentation](http://objectpath.org/reference.html).

#### Transparent projections
Before we have only encountered two types of prototypes: directory and file. There are some cases when we need to flatten 
complex nested resource structure into more manageable structure. In this case we have prototype type called: transparent prototype.

Let`s create projection to illustrate transparent prototypes capabilities. Firstly we remove old filesystem projection example 
to free up mount point (or we can create new mount point and project second example into it):

```
$ ./projections_cli.py rm -n fs_example_1
```

Now we create example flat projection:

```
$ ./projections_cli.py project -n fs_example_2 -m mount -p examples/fs_example_2_config.yaml -d fs_driver
```

After projection initialization we will see this structure using "ls -R mount":

```
mount/:
test_dir

mount/test_dir:
sample_1_bed_file.bed  sample_2_bed_file.bed  sample_3_bed_file.bed  sample_4_bed_file.bed  sample_5_bed_file.bed
sample_1_vcf_file.vcf  sample_2_vcf_file.vcf  sample_3_vcf_file.vcf  sample_4_vcf_file.vcf  sample_5_vcf_file.vcf
```

In projection configuration file on path "examples/fs_example_2_config.yaml" we can see that "sample" prototype 
definition is not so different from "sample" prototype of previous example. Apart from rerun directory prototype, 
excluded for brevity only "type" field is different. Children nodes too have change only a bit, but this is 
important change in field "name" of "fasta" prototype, for example:

```
name: " join([$.environment.name, $.name], '_') "
```

Here, projection name is joined using "sample" directory name and "fasta" file name. This step is essential to avoid 
names collision on one level of projection, **there must be no equal paths in one projection**.

#### Binding metadata

Metadata - data relations in projections meta filesystem are defined as links between projection nodes. Each link has 
it`s head node (node which will be annotated with metadata) and tail node (node which is metadata for head node).  