# A little bit on Clustering

#### Skills: Familiarity with sklearn package and an understanding of why clustering is useful in bioinformatics

#### Preliminary note:

This lesson will force you to look through documentation and fill in the blanks on some of the functions. I know I could have given you the full code, but I think that more is learned when you are forced to interact with the code to some extent. Note that in Python syntax, you will have to use the dot (period) to indicate when you are going from a general to specific scope. For example, there is a package called numpy, with a module called random, with a function called seed. The seed function is called by typing out ```np.random.seed()```, telling python where exactly it needs to look. 

**WE ARE USING PYTHON3, NOT PYTHON**

### General Factoids

**Clustering:** The process of partitioning a dataset into groups based on shared unkown characteristics. In other words, we are looking for patterns in data where we do not necessarily know the patterns to look for ahead of time. 

**In:** N points, usually in the form of a matrix(each row is a point). 

**In:** A distance function to tell us how similar two points are. The most intuitive distance function is [euclidean distance](http://rosalind.info/glossary/euclidean-distance/). The type of distance measure you use can vary depending on the type of data you are using. One specific example of how fine tuned distance metrics can be is a distance metric that weights substitutions in DNA strings heavier than indels because indels are more common sequencing errors. This means that such a metric would not separate two DNA strands due to sequencing error. 

**Out:** K groups, each containing points which are similar to each other in some way. Some algorithms cluster for a preset number K, others figure it out as they go along. 

**How is DNA a point?** Discussing points in space is nice, but a little abstract. The way we translate DNA to points is by making a string of DNA into a kmer vector. A kmer is a string of length k, and there are 4^k possible kmers in any DNA strand. Usually we think of kmers as ordered lexicographically like this: AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA and so on. One way to represent a DNA strand is to create an array of zeros of length 4^k and increment by one for each time that kmer appears in the DNA strand. 

*Problem:* Find the kmer vector of "AGTTTCAT" for k=2. 

### Example Bioinformatics Applications: 

1. Finding what genes are up and down regulated under certain conditions. Imagine you have a matrix, where each point is a set of gene expression recorded for a variety of conditions. If two points are close to each other, that means they had similar expression levels throughout those conditions. 

2. Discerning different species present in a sample of unkown contents. This can be done with an algorithm that does not have a preditermined amount of clusters. An example application is taking HIV sequences from a patient, clustering them, and filtering clusters under a certain size to find the sequences of prevalent strains within the patient.

3. Finding evolutionary relationships between samples using hierarchical clustering. The earlier on two centers were combined, the closer their corresponding points are from an evolutionary perspective. 

## Clustering Methods

In this class, I want to demonstrate the difference between different clustering algorithms. To start with, we need data that will make these differences very obvious. 

First off, copy the file ```/home/ubuntu/clustering_lesson/clust.py``` into your own directory. 

Start by creating some blobs using numpy's random normal generator. First, create a seed for the random number generator with ```np.random.seed()```. Then, use the following code as a guide for making a blob: ```clust1 = np.random.normal(5, 2, (1000,2))```. This will make a blob centered around 5 with a stdev of 2 in the form of a 1000\*2 array.

Make a few blobs (and save them in variables `clust1`, `clust2`, etc.) with centers between 0 and 20 and varied stdevs (keeping the array dimensions the same). 

Now, use ```set1=np.concatenate()``` to bring the blobs into one structure. (**Hint**: Don't know how to use this (or any other function below)? Try googling to find the official documentation!)

The figure we can use to contrast the blobs is some circular point layouts. Make two concentric circles with ```set2=datasets.make_circles(n_samples=1000, factor=.5, noise=.05)[0]```

Finally, use the provided ```cluster_plots()``` function with the two sets you have generated to generate some pretty plots. **Note:** this will save your plots to a pdf file in your current directory, the name of which can be changed (look at the argument called *name*!) You'll have to scp the file over to your computer to view it. 

<details>
  <summary>Forgot how to scp? click here</summary>
  
```
scp username@ec2-18-218-72-125.us-east-2.compute.amazonaws.com:/home/username/path /localpath/
Do not forget to replace the paths and username!
```

</details></br>

### 1. K-means

&nbsp;&nbsp; I. Select K centers. This can be done a variety of ways, the easiest of which is to select a random point, find the farthest point from that and select it, find the furthest point from the previous and so on. 
  
&nbsp;&nbsp; II. Assign each point to the center nearest to it. 
  
&nbsp;&nbsp; III. For each center, take all the points attached to it and take their average to create a new center for that cluster.
  
&nbsp;&nbsp; IV. Reassign all points to their nearest center and continue reassigning + averaging until the iteration when nothing changes
  
  Code to apply K-means clustering to set1 is below, don't forget to replace k with the number of blops you made above. Do the same for set2, but now set k=2.   There are two circles, right? It seems like each circle might share a characteristics within itself, so it would be nice if the computer could cluster the circles into 2 partitions. 
  ```
  kmeans_dataset1 = cluster.KMeans(n_clusters=k).fit_predict(set1)
  cluster_plots(set1, set2, kmeans_dataset1, kmeans_dataset2)
  ```
  
### 2. Agglomerative Hierarchical 

&nbsp;&nbsp; I. Start with every point being its very own cluster center. 
   
&nbsp;&nbsp; II. Find the two centers which are closest to each other and combine them into one cluster. The new center is the average of the two previous ones. 
   
&nbsp;&nbsp; III. Continue combining the closest centers until all points are under a single cluster. 
   
   Let's look at the blobs and circles we made again. Look at the [documentation page for agglomerative clustering](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html) and figure out the proper calls for the two datasets (syntax is very similar to what you did for k-means).
   
   Now, run `cluster_plots()` to graph and `scp` to view the results. Do you notice a difference in quality of the clustering of the left and right graphs?
   
   Unfortunately, this still does not solve our circlular cluster issue. For that, we can ask the sklearn package to build a graph out of our circles which restricts the amount of nearest neighbors a point can cluster with. 
   
   ```
   #this line is missing a parameter - how would you set the number of nearest neighbors to 6?
   connect = kneighbors_graph(dataset2, include_self=False)
   
   #in this line, you need to set linkage to complete, number of clusters to 2, and set connectivity equal to the graph 
   #on the previous line
   hc_dataset2_connectivity = cluster.AgglomerativeClustering().fit_predict(dataset2)
   ```
   
   Use ```cluster_plots()``` to graph the plot resulting from the regular AgglomerativeClustering beside the graph that plots the connected AgglomerativeClustering. As you may have guessed, this fixes the circle problem. 
   
### 3. Soft Clustering
  
  I. Select K centers randomly (or slightly less randomly, like [farthest first traversal](https://en.wikipedia.org/wiki/Farthest-first_traversal). 
  
  II. Expectation step: Assign a responsibility (a likelyhood) for each point in respect to each cluster. Basically, the closer a point is to a center the higher the likelyhood of that point. This can be calculated in a variety of ways, but the main idea is that it is some function relating the distance from a point to a center to the distances between all points and that center to see if it is much closer or further than other points
  
  III. Maximization step: compute new cluster centers by using a **weighted** average of the points based on the likelyhoods calculated in step II. 

  Let's go ahead and see what will happen with our blobs and circles under this clustering algorithm:
  
  ```
  #I think you can figure out the arguments we need for expectation maximization
  em_set1=mixture.GaussianMixture().fit_predict
  ```
  
  Use these new clusters in our ```cluster_plots()```


## Dirichlet Process Means Clustering

I provided a starting file at ```/home/ubuntu/clustering_lesson/dmp.py``` to get the ball rolling 

Now that you have learned some terminology so you know how to Google information about clusters, let's practice Python by making our own clustering algorithm. 

The dirichlet process is just a method of clustering without knowing the number of clusters ahead of time. Something like this algorithm may be used for example 2 in the list of clustering related bioinformatics problems. Here are the steps:
 
&nbsp;&nbsp;  I. Add the first point in your data as the first center
  
&nbsp;&nbsp;  II. Iterate through all of the points in your data, calculating the distance between the current point and all existing centers. 
  
&nbsp;&nbsp;&nbsp;&nbsp; A. If the point falls within a certain threshold distance of the center, add that point to that center's cluster
  
&nbsp;&nbsp;&nbsp;&nbsp; B. If the point does not fall within a certain threshold distance of the center, add that point to the list of centers
    
&nbsp;&nbsp; III. Once you have gone through the Dirichlet Process, you do the means part. Redifine the centers of each cluster to be the average of all points in the cluster, like you would in most standard clustering algorithms. 
  
&nbsp;&nbsp; IV. Repeat steps II and III either until an iteration provides no change in the assignment of points to centers, or a maximum number of cycles is reached. 

**HINT:** I recommend using arrays to store the cluster information! You will need an array that is the same length as the number of points(index 1 is point 1 and so on). Then, you will assign a cluster number to each index as you go along (keep track of the number of clusters//increment when you create a new one). Also, you need an array of the actual centers. 
  
