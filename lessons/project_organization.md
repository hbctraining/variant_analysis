# Project Organization

## Learning Objectives
- Configure a workspace on the `scratch` drive
- Organize dataset for analysis
- Differentiate between using `home` and `scratch` drives

## Logging into O2

For this workshop we will be using training accounts to log in. These have been created for us by the [HMS Research Computing (HMS-RC) team](https://it.hms.harvard.edu/our-services/research-computing), they are the folks that manage the O2 cluster. We will be providing each of you with your own training account associated with a password for the duration of this workshop. Your training account and password can be found [here]().

> If you are interested in getting your own personal account on O2, please follow the instructions provided here after this workshop.

Let's get started with the hands-on component by typing in the following command to log in to O2:

```
ssh username@o2.hms.harvard.edu
```

You will receive a prompt for your password, and you should type in your associated password; ***note that the cursor will not move as you type in your password***.

A warning might pop up the first time you try to connect to a remote machine, type "Yes" or "Y".

Once logged in, you should see the O2 icon, some news, and the command prompt, e.g. [rc_training10@login01 ~]$.

> Note 1: ssh stands for secure shell. All of the information (like your password) going between your computer and the O2 login computer is encrypted when using ssh.

## `home` or `scratch`

Within the [Introduction to the Command-line course](https://hbctraining.github.io/Intro-to-shell-flipped/schedule/), we introduced your `home` and `scratch` workspaces, but let's give a brief recap of each and their purpose:

`home`
 - Limited to 100GiB of storage per user
 - Backed-up daily up to 14 days, weekly up to 60 days
 - Can you found at `/home/$USER/`
 - Great place for storing valuable file/scripts

`scratch`
 - Up to 10TiB or 1 million files/directories per user
 - No back-ups
 - Purged from the system if not accessed for 30 days
 - Great for storing reproducible, intermediate files

Due to the limited storage space on `home`, we are going to take advantage of `scratch` to hold some of our intermediate files for this workshop. This is a very common use of the `scratch` space as many analyses will have large intermediate files, which would otherwise fill up our `home` directories.

### Creating a `scratch` space

While on the login node, we will create our space on `scratch`. In order to create a space on `scratch`, we will need to run a script provided by the HMS Research Computing team:

```
$ sh /n/cluster/bin/scratch3_create.sh
```

> Note: You *MUST* be on a login node in order to create a space on `scratch`.

It will prompt you with the following:

```
Do you want to create a scratch3 directory under /n/scratch3/users? [y/N]> 
```

To this you will respond `y`, then hit `Enter/Return`.


Next, it will prompt you with:

```
By typing 'YES' I will comply with HMS RC guidelines for using Scratch3.
I also confirm that I understand that files in my scratch directory WILL NOT BE BACKED UP IN ANY WAY.
I also understand that THIRTY DAYS after I last access a given file or directory in my scratch directory,
it will be DELETED with NO POSSIBILITY of retrieval.

I understand HMS RC guidelines for using Scratch3: 
```

Type `YES`, then hit `Enter/Return`.

It should return:

```
Your scratch3 directory was created at /n/scratch3/users/r/rc_trainingXX.
This has a limit of 10TiB of storage and 1 million files.
You can check your scratch3 quota using the scratch3_quota.sh command.
```

Now that we have created our `scratch` space, you will need to start an interactive session. A login node's primary function is to enable users to log in to a cluster, it is not meant to be used for any actual work/computing. Since we will be doing some work, let's get on to a compute node:

```
$ srun --pty -p interactive -t 0-3:00 --mem 1G  /bin/bash
````

Make sure that your command prompt is now preceded by a character string that contains the word `compute`.

## Organizing our project

The first thing we should do is make sure that we are at our `home` directory by using:

```
cd ~
pwd
```

This should let us know that we are located at: `/home/rc_trainingXX`

Once we have established our location, let's begin by making the directory that we will be using to hold our work in the `home` workspace:

```
mkdir variant_calling
```

Next, we are going to move into this newly created directory and create some additional directories to hold our work:

```
cd variant_calling
mkdir scripts results figures raw_data reports
```

Let's move into our `raw_data` directory and copy the FASTQ data that we are going to use for our analysis:

```
cd raw_data
cp /n/groups/hbctraining/variant_calling/raw_data/*.fq.gz .
```

This may take up to a minute as it has a lot of data to copy. Now that we have created the directories that we are going to use in our `home` space. Let's move to our `scratch` space so that we can set up the directories that we are going to there to hold our intermediate files:

```
cd /n/scratch3/users/${USER:0:1}/$USER
```

A few quick notes about this previous `cd` command:

1. You can see that we use the `$USER` variable twice in this `cd` command. `$USER` is an environment variable that is your username. You could check it by using `echo`:

```
echo $USER
```

2. Perhaps the more puzzling part of the command is the use of `${USER:0:1}` in the path. The `users` on the scratch space are organized into subdirectories starting with the first letter of their username. `${USER:0:1}` will return the first letter of the username. You can test it using `echo` as well:

```
echo ${USER:0:1}
```

By using `${USER:0:1}`, we have seen how you can produce substrings in `bash`. The syntax is:

```
${VARIABLE:START:LENGTH}
```

Where, 
  - `VARIABLE` is the variable that you would like to create a substring of
  - `START` is the 0-based indexing of the position to start the substring
  - `LENGTH` is the number of characters following the starting position to include in the substring

Of course, it you knew the first letter of your username you could use that instead of `${USER:0:1}`. We choose to use `${USER:0:1}` here and throughout the rest of the lesson in order to hopefully make this code more easily reuseable for your own research.

Now that we are located in our `scratch` space, let's go ahead and create a directory for our analysis:

```
mkdir variant_calling
```

Now, move inside this newly created directory and create the following directories to hold our intermediate files:

```
cd variant_calling
mkdir alignments vcf_files
```

Excellent! Our project space has now been organized and is ready for the analysis to begin!

[Next Lesson >>>](file_formats.md)

[Back to Schedule](../schedule/README.md)


