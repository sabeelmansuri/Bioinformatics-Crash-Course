# The Linux Terminal 

#### Skill: UNIX/Command Line + RNA-seq quantification

The previous doc we went over included some information on alignment. This lesson will teach command line skills which will allow us to use alignment software. There are an endless number of commands, each with a ridiculous amount of options, so **do NOT attempt to memorize on the first try**. Use the commands listed below as a reference to look at. Actually learning (or memorizing) the commands comes from repeated use of the terminal. Same goes for the tool installation process we will go through. This course is supposed to be all about these tools, so we will have much more practice with handling them. Think of this as an introduction, a reference sheet, and generalized example.

#### About the RNA-seq Tasks

As a way of practicing command line skills, this lesson will also go through the motions of a very common problem in bioinformatics: quantifying gene expression in RNA-seq data. Today we will be looking at a specific type of RNA, called lnc RNA. Long Non coding(lnc) RNA, strands of 200 or more non coding nucleotides are implicated in a range of gene regulation roles, from epigenetic X inactivation to transcriptional control. The irregular expression of lnc RNA is involved in the parthenogenesis of just about every cancer, including gastric cancer.

Due to data storage restrictions, I have subsampled the RNA-seq data for you and placed it in ```/srv/lesson2/sub_Gastric_Ctr``` and ```/srv/lesson2/sub_Gastric_Affect``` (sub stands for subsampled). The RNA-seq data was obtained through a paired-end illumina sequencing protocol. This means that there is a forward and a reverse read for every fragment which was sequenced. According to illumina, this makes for better alignment to a reference genome. You will find there is a file called ```sub_SRR2073159_1.fastq``` and ```sub_SRR2073159_2.fastq```, which are the forward and reverse reads of the same fragments. Software made for illumina data processing will always have a paired-end option, so pay attention to what kind of data you have! If you want to look at where I got the data for today's lesson, [go here](http://cancerpreventionresearch.aacrjournals.org/content/9/3/253) and take a look at the "Next-generation sequencing analyses" section. It will give you an accession number corresponding to the experiment's data. Accession numbers will be a story for a future lesson(downloading data takes horridly long so I'm going to have to think of a way to get around this while teaching the subject). 


## But first... EC2

I set up accounts for everyone based on the usernames from last time. Click [here](https://docs.google.com/spreadsheets/d/1M4S22RieI7GnJqGJZo_4flSU3FzP7ypCqrNSjZ-rf9w/edit?usp=sharing) for the list of usernames.

***Secure Shell(ssh):*** a protocol which creates a secure channel for two computers to communicate even over an unsecured network. This is how we will connect to EC2. 

**Windows:** Open putty, paste ```my_username@ec2-18-188-27-13.us-east-2.compute.amazonaws.com``` into the Host Name section and select port 22 and SSH on that same page. Type a name under Saved Sessions and click the Save icon on the right. Now, press Open and type your username and password when prompted. 

**Mac or Linux:** Right click anywhere and click open terminal. You should see a prompt that looks something like  ```mchernys@mchernys-ThinkPad-T430:~/Desktop$```. Next, copy paste this command into the terminal and press enter ```ssh my_username@ec2-18-188-27-13.us-east-2.compute.amazonaws.com```. Note: use Ctrl-shift-V to paste into terminal. Please replace "your-username" with your actual username. Save the command you used somewhere so you can copy paste it in the future. 

**Some things to keep in mind:** I have some promotional credits for EC2, but they are not infinite. Please help me not run over budget by following my instructions and asking if you are unsure of what you are doing. Processing power is plentiful and the rate charged is constant, but storage can add up if everyone uploads large files. Please only upload what I ask you to upload. A few megabytes here and there is fine but please do not go uploading several gigabytes at a time. In addition, the actual process of uploading costs no money but downloading from EC2 does. Again, just download when I ask you to. Thanks!!

---

## Passwords!

Right now, everyone's passowrd is ubic2018 but you probably have your own password in mind for your account. Type ```passwd my_username``` and follow the prompts to set your own password. 

---

Discover your identity. Type `whoami` into the window that just opened up and hit `enter`. And just like that you're talking
with your computer, you bioinformatician, you.

## Why should I learn this?

Quoting Mark's PI: "ja, true, you need that shit."

Other reasons: It is the fastest way to deal with XXL sized files, much of the software in bioinformatics has no Graphical User Interface (GUI), and any programmer should have a healthy relationship with their system.  

## How do I see what a command does?

Anytime you need a refresh on what a command does, type the command line with the --help option like so: ```ls --help```. If that does not work, try ```man ls```. I will go over why different commands have different help syntax in a bit. 

## Navigation, manipulation, and permission

In order to get started, we need to be able to do the same thing we do in a file explorer in the command line. You may find it inconvenient at first, but with time these commands become faster and more versatile than the file explorer's interface. 

The forward slashes in a terminal console represent directories, with the home directory being a ```~```. Your default folder on EC2 is your use folder, which is ```~/username```. This means the folder named after your username is a subfolder of the home folder, which is represented by ```~```. 

```screen``` A way to have multiple terminal instance from a single screen. Start screen by typing ```screen -a```, open new screens with ```ctrl a c```, move between screens with ```ctrl a n```. If you already have screeens working in the background, you can access them with ```screen -dr``` (so you don't have to start with a blank terminal every time). 

*Note: the ```ctrl d``` command for closing terminals will close one screen window at a time, as one would (hopefully) expect it to do* 

```cd```(change directory) Type cd followed by the directory's path to navigate a terminal to that directory. ```.``` is current directory and ```..``` is the parent of the current directory. 

```ls```(list files) prints out the contents of a directory. There are tons of options for this command - my favorite is ```ls -lah``` , since it prints the directory contents in list format(```-l```), includes hidden files/folders(```-a```), and makes the storage sizes more readable for humans(```-h```). 


```mkdir```(make directory) Creates a directory with the same name as the argument you give it. 

---

### TODO: Make a Software Folder

Navigate your terminal to your home directory (the directory named after your UCSD username) using ```cd```. Type ```mkdir software``` and press enter. Type ```ls``` to see the changes you have made. The reason for a software folder is to... keep your software in it. 

*Note: Usually, you would place executables in the /bin system folder, but you are not the admin so you cannot access that folder :( . This is often the case when you ssh into a system, so get used to having a dedicated software folder.*

---

```cp```(copy) copies the file in the first argument to the directory in the second argument. ```cp file1.txt file2.txt``` makes a copy of file1.txt called file2.txt in the same location. ```cp file1.txt ..``` places a copy of file1.txt (called file1.txt) into the parent directory. 

```rm```(remove) deletes a file. Careful with this one, there is no convenient command to remove files from the trash on linux. Once deleted, gone forever. ```rm file.txt```

```mv```(move) like copy, but the original file disappears. 

```wc```(word count) counts things like lines, words, and characters. ```wc -l file.txt``` prints the number of lines in file.txt. 

```chmod``` In order to execute files, you need permission to do so. When looking at the output of ```ls -lah``` , you will see something on the order of ```-rwxrw-r--```. This indicates that the owner has read, write, and execute permissions. The next three characters are group permissions, and the last three are permissions for everyone else. Change permissions with ```chmod XXX filename``` where each X is a number 1 through 7 (first for owner, second for group, third for everyone else).

***Permission	rwx	Binary***

7	read, write and execute	rwx 111

6	read and write	rw-	110

5	read and execute	r-x	101

4	read only	r-- 100

3	write and execute	-wx 011

2	write only	-w- 010

1	execute only	--x	001

0	none	---	000

## Tying It All Together 

### Piping 

Piping is stringing multiple commands together to perform some larger task. You can string commands together using the `|` symbol like so:

```
cat file.txt | grep "hello"
```

This will first perform the `cat` command: it will attempt to print everything inside file.txt to the terminal, but the pipe symbol `|` will stop it. Instead, everything from file.txt will get pushed through the `grep` command, which will print out any line that contains "hello". The result is only the lines containing "hello" in file.txt will be printed out.

## Getting yer feet wet

Here are a quick batch of tasks/exercises using the command line (and different commands) in roughly increasing order of difficulty:

##### (0.) Copy-paste the following two commands into your terminal in order and hit enter.
`cd ~`  
`cp -r ../smansuri/parent/ ~`

##### 1. Enter the new directory (hint: use `ls` and `cd`)

##### 2. Name all of the files (not other directories!) inside this directory. How many are there? (hint: use `ls`)

##### 3. How many lines are in file1.txt? (hint: use `wc`)

##### 5. How many lines inside file2.txt contain the characters (not the word) "no"? This means "anNOunce" will count too. (hint: use a pipe of `grep` and `wc`)

##### 6. (Challenge 1) Execute the "instructions" file. Follow the instructions. (hint: `./instructions`)

##### 7. (Challenge 2) Attempt to delete the directory you downloaded in step 0. This directory is called "parent".


## Downloading

```scp```(secure copy) is a command used to copy files from one machine to another. The first argument is the source location, while the second argument is the destination. ```scp file.txt my_username@dns_address.com:/home/my_username/docs```

```curl``` Will download stuff for you. The most simple and relevant combination of options is ```curl -L https://examplelink.com -outdir .``` which will download from https://examplelink.com into the current directory (indicated by the dot). 

```apt-get```Handles packages from the apt library for Debian based systems. However, this installs packages system-wide so you are not going to be able to use it on EC2. The mac equivalent is homebrew. ```sudo apt-get install google-chrome-stable``` will install chrome. 

---

### TODO: Get FastQC and Kallisto

Quality of genetic information is important! FastQC is the gold standard for quality control in the bioinformatics field. Google "download FastQC" and find the instructions. Make sure to select the correct package for a linux system, then go think about how you would get that package onto EC2. There are two main ways to do this. 

*Hint: the last three commands mentioned contain both of the two ways you can get FastQC onto EC2*

Next, we need Kallisto, which describes itself as "a program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads." Google "download kallisto" and find the appropriate file (it should be a .tar.gz). Use the same method you used for FastQC to transfer it to EC2(or challenge yourself to find the second way).

---

## Unpackaging

Much of the data people want to download is large, but they want it fast. That's why things like .zip, .tar, .gz and such exist. Those are the file extensions of compressed data. In order to make software work, it must be unpackaged.

```unzip``` Is exactly what it sounds like. This command unzips .zip file types. 

```tar```(tape archive) Is the command linux uses to package and unpackage stuff. This command has an incomprehensible amount of confusing options, so let me just copy paste the ones you should care about. ```tar -xvf file.tar.gz -C .``` unpacks a .tar.gz file into the current directory and ```tar -xvf -C .``` unpacks a .tar file into the current directory. The -C option indicates the files' destination.

---

### TODO: Unpackage your FastQC Kallisto

Now you have your fastqc*.zip in your software folder. In order to use it, you're going to have to unpackage it. Look a couple lines up to figure out how. 

The same goes for your kallisto .tar.gz. This is a filetype you will run into often when dealing with linux, since it is the default way linux compresses a folder. If you look in the Unpackaging section of this document, you will find instructions on how to open up this strange creature. 

---

## Compilation

What is compilation? It is the conversion of one programming language into another. Typically, it is a conversion of what's known as a high-level language (C, Java, Python, etc) to a low level language (binary, assembly). CPUs understand only very very very basic logic, so a super smart program called a compiler has to convert your convoluted and messy code into the simple delicious porridge that the CPU can eat(execute).

**shell scripts and python files** do not need compilation.

**java** compiles by ```javac filename.java```

**C** ```gcc -c filename.c``` to compile and assemble. 

Note: Much of the time, software you download online is already in binary form so there is no need to compile. This is not always the case!

## Execution

**/bin**(binaries) contains your executable files and shells. The computer has a list of folders it searches through to find executables when you type a command and the /bin directory is one of them. When you download software, you should place the executable file or a symbolic link into the /bin directory.

**shell script** a simple ```./executable``` will suffice to execute a script. 

**python** will automatically compile for you before executing with the command ```python filename.py```. 

**java** can be executed with ```java compiledfilename```

**C** is executed like a shell script ```./out```

---

### TODO Actually Quantify Stuff

Okay, we now how to execute now. The FastQC folder you have now contains an executable called fastqc. Go into the folder containing the executable, type ```./fastqc --help``` to see usage instructions. Your task is simply to run the fastqc on each of the files sitting in the ```/srv/lesson2/sub_Gastric_Ctr``` and ```/srv/lesson2/sub_Gastric_Affect```directories. FastQC will produce some .html files, which I will just show on the large screen and explain in the interest of saving time. 

Kallisto is a little more complicated. Our first step is to build a kallisto index file, which will assign an index to each RNA transcript we will quantify and optimizes the quantifying procedure in general. Where did we get these RNA transcripts you ask? Good question! I googled "mouse lncRNA database" and eventually happened upon a [website](https://www.gencodegenes.org/) where there was a fantastically easy to download .fasta file full of lnc RNA sequences from mice. Anyways, go back to the kallisto website and look at the instructions on how to index a file. The list of target sequences I got from gencode is at ```/srv/lesson2/gencode.vM17.lncRNA_transcripts.fasta``` and you might want to specify the name for the destination file. 

Next, you will need to look at the ```kallisto quant``` command. Specify an output folder, the index you made just now, set --bootstrap-samples=100, and finally include the forward and reverse paired end files at the end. You will run the quant command 4 times for the 2 control and 2 affected files. You might want to create 4 separate output folders in order to keep everything organized. 

The last part will be a bit of a walkthrough, since it is kind of complicated (I don't understand every option either, don't worry about it). Open up R  by simply typing R into the terminal and pressing enter. Below is the template for what you need to execute in R in order to compare transcription levels. Note the parts where you need to enter a path and replace them with your own filepaths


```
library("sleuth")
#paths to Kallisto outputs from MPNST study
sample_id = c("/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Ctr_Rep1", "/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Ctr_Rep2", "/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Affect_Rep1", "/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Affect_Rep2")

kallisto_dirs = file.path(sample_id)
s2c = read.table(file.path("/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/sleuth_Gastric_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c = dplyr::mutate(s2c, path = kallisto_dirs)

so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)


write.table(sleuth_significant, "/home/my_username/sleuth_output/", sep="\t", quote=FALSE)
gg=plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log", offset = 1)
ggsave("/home/my_username/plots/", plot=gg)
```
