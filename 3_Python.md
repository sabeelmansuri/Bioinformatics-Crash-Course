# Intro to Python
#### Skill: Python

## Preface

### If this is your first time here, please check in with Mark to get an AWS account. Then, a board member can help you connect (SSH) to the cluster

### Please fill out this sign-in sheet: https://goo.gl/forms/YoNHEmLpsYnFl0A62

### The cluster is live on: [username]@ec2-3-16-76-94.us-east-2.compute.amazonaws.com

## The Big Picture

So far, we've been working with bioinformatic data only on the command line. However, you'll often want to do more than the command line easily allows. Today, we're going to be learning the basics of one of the most popular languages for doing that: Python.

(Side note: If you know a bit about Python, you'll know that there are two commonly used versions: python2 and python3. **They are VERY similar.** You're going to run into programs written for/in either version, so being familiar with both is important. For this lesson, however, we're going to be using python2 because... well... we are.)

## Getting Started

Create a new directory (remember: mkdir) wherever you are saving your work. Call it "pydir". Enter the directory.

Python should already be installed on your workstations. Let's make sure: Type the following on the command line:
```shell
python
```
You should see something that looks like this:
```shell
Python 2.7.15rc1 (default, Nov 12 2018, 14:31:15)
[GCC 7.3.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>>
```
Type exit(), and move on.

## How to Python

### Hello World!
We're going to start by taking a quick look at how Python code is written. Create a file called Hello.py. The ".py" extension signals that you will be writing in python. Inside the file put:
```python
print "Hello World"
```

Great! Now save + close the file, and run your newly-written program by typing the following on the command line:
```shell
python Hello.py
```

Because the "print" statement in Python outputs whatever follows it to the command line, you'll see your program print "Hello World". That was pretty trivial... let's try something more interesting.

### Basic syntax

You'll need to obtain a file I've written and add it to your "pydir" directory. Using the `cp` command, and knowing that the file is currently located in `~/../smansuri/syntax.py`, copy the file into your "pydir" directory.

Now, open the syntax.py file and take a look inside.

Make sure you understand what's happening. Follow the comments closely, and ask one of us if you have any questions.

Now, run the code the same way we ran "Hello.py". You should see "Bioinformatics is Cool". Can you edit line 14 to make the program print: "is Bioinformatics Cool"?

##### vim-hint: navigate to a specific line in vim by going into command mode and typing ```:14```

At this point, you can remove the "#" from the start of the last line (this is called uncommenting). Your Python senses should tell you that this line will now print out myString. Take a look at how myString is defined above and take a guess about what should be printed when you run the program. Once you're ready, run the program.

What if you only wanted to print *part* of your string, not the whole thing? Remember that a string is like a list (with the first character at index 0). So, what if we wanted to print just "Hello"? We can use specify a range of indexes to print from like so:

``` python
print myString[0:5]
```

Note that the first number is the position of the first character printed (0 = 'H'), while the second number is **PAST** the last character printed (5 = ' ', but we only print up to index 4).

## Indentation 

Indentation in Python **matters**. Try adding a second print statement your "Hello.py" file so it looks like this:
``` python
print "Hello World"
    print "Indented line"
```

Now, try to run "Hello.py". Python will complain that there's a problem with your indentation (there was no need to indent, but you did anyways). You'll learn more about when to indent in the next section. Speaking of which, it's time for some bioinformatics.

## Loop-D-Loop

Say I have a fasta file that has some number of reads (make sure you remember fasta format or you'll have trouble with this exercise!). I want to write a Python program that ONLY outputs the header lines (the ones that start with ">"). How can I do it? Think about how you might start before continuing.

```
-check every line-
  -if it starts with a ">"-
    -print the line-
```
This is one simple **representation** of how you could achieve this task. The **implementation** in Python, as we shall see, uses a loop. Which one of the three pseudocode lines above suggests we will need a loop?

First, however, let's get a sample fasta file. Use the following command to download it straight in your working directory:

```shell
cp /srv/Winter19/gencode.vM17.lncRNA_transcripts_subsampled.fasta test.fasta
```

Check the contents of your directory. You should see a file called "test.fasta". 

Let's make a Python program that reads from this file. Create a new file called "Loop.py" and add this as the first line:
```python
file = open("test.fasta", "r")
```

This will open the file (test.fasta) for "r"eading, and give you access to test.fasta in a variable called "file". Now let's use a loop to look at every line in the file:
```python
file = open("test.fasta", "r")

for line in file:
```

That last line is the syntax for starting a for loop in Python. Next, looking back at our pseudocode, we see that we need to check if a line starts with a ">". Luckily, lines in a file are stored as strings! Remember that strings are indexed, meaning individual characters from them can be extracted using brackets (you might remember we did something similar above with!). 
```python
file = open("test.fasta", "r")

for line in file:
  if line[0] == ">":
```

Pay close attention to the indentation here. **You can think of everything that's indented after the `for` as being "inside" of the for loop** (it looks a bit like that too!). In our example, that means the `if` code is executed for every line in the file.

Finally, we just specify that we want the line printed if the line does start with ">". We indent the next line so Python knows it's part of the "if" statement, and...

```python
file = open("test.fasta", "r")

for line in file:
  if line[0] == ">":
    print line
```

Run Loop.py and see what happens. Voila.

## Packages

An important feature of Python is its packages. Think of packages as Python programs someone else has written that you can borrow and use for your own purposes.

For example, you may want to take the first sequence of "test.fasta" and use NCBI BLAST (learn more about the BLAST database [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) to determine what organism the sequence may have come from. You wouldn't want to write your own code to connect to the database... instead you can use a nifty package to do that for you. Let's try:

First, extract the first two lines of "test.fasta" into a new file, "small.fasta":

```shell
head -n2 test.fasta >> small.fasta
```

Now, install the correct package. The one you want to use is called [BioPython](https://biopython.org/). Use the Python program installer (called pip) to install the package to your user:

```shell
pip install --user biopython
```

Great! BioPython is now available on your account. Let's use it. Create a new file "Blast.py", and add the following code:

```python
from Bio.Blast import NCBIWWW

fasta_string = open("small.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string, format_type='Text', hitlist_size=1)
print result_handle.read()
```

Let's walk through what this does. The first line takes the programs we want for NCBI from BioPython and prepares them to be used. The second line reads in the "small.fasta" file. The third line is the most important: it takes the genetic data, connects to the NCBI BLAST database, searches for matches, and then returns the result from the database. Finally, the last line prints that result.

You can try running "Blast.py" now if you'd like, but I'd recommend coming back to this after completing the next 3 exercises.

## Your turn!

### Warm-up
You're going to print your name using the alphabet. Declare a variable that holds a string containing
all 26 letters of the alphabet, and one space character (just use the declaration provided below for simplicity--the last character is a space!). Then, on the next line use one print statement that accesses letters in the alphabet variable to print out your name.

```python
"ABCDEFGHIJKLMNOPQRSTUVWXYZ "
```

We'll now shift our focus onto applying what we've learned to bioinformatics. 
Various implementations may work. Pick any that you like.

### Challenge 1

Create a file that takes in the DNA data (A/T/C/G) from "test.fasta" and prints out the file as RNA data (A/U/C/G). (Note: This is not how real transcription works. Just use this simple replacement as an example).

**Make sure you don't attempt to transcribe the header lines!**


### Challenge 2

Create a file that takes in "test.fasta" and prints out the GC-content of all of the data. Hint: the GC-content is the percentage of G's + C's in some genetic data. If you wish, you can read more about it [here](https://en.wikipedia.org/wiki/GC-content). Don't worry about whitespace stripping (especially if you don't know what that is).

**Hint: If you keep getting 0, try explicitly converting either the GC count or the string length to a decimal using float()**

## Credits

Exercises were adapted from [Rosalind](http://rosalind.info)
