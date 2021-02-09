# Nucleus Decomposition of Probabilistic Networks
This repository contains efficient implementations for computing nucleus decomposition in probabilistic graphs. The link to our dataset can be bound in https://www.dropbox.com/sh/cunq076ocfu8m05/AACU5Wcz9p_vp5PJ-0RsBbdja?dl=0

# Local Nucleus Decomposition
1) **PNucl_DP_Efficient.java:** This is an implementation of local nucleus decomposition of probabilistic graphs which uses dynamic programming for computing and updating the support of triangles (DP).

2) **PNucl_CLT_BinoJava_DP_Efficient_RegularizedBeta.java:** This is an implementation of local nucleus decomposition of probabilistic graphs which uses statistical approximations for computing and updating the support of triangles (AP).

# Weakly-global Nucleus Decomposition  
**Weakly_Global_onSmallGraph_ForExperiments.java:** This is an implementation of weakly-global nucleus decomposition of probabilistic graphs which gives approximate solutions based on search space pruning combined with Monte-Carlo sampling.

# Global Nucleus Decomposition
**Golbal_nuclei_Finding_onSmallGraph_ForExperiments.java:** This is an implementation of global nucleus decomposition of probabilistic graphs which gives approximate solutions based on search space pruning combined with Monte-Carlo sampling.

# Input 
The graph datasets used for our algorithms should be in WebGraph format with edges being assigned probabilities. We refore to this type of Webgraph as ArcLabelled Webgraph. 

There are three files in this format:

newTest.w.labeloffsets<br/>
newTest.w.labels<br/>
newTest.w.properties<br/>

see newTest example in the main directory. 

There are many available datasets in http://law.di.unimi.it/datasets.php which can be converted to an ArcLabelled Webgraph. These datasets are unweighted and directed graphs, which are in Webgraph format.

Let us see for an example dataset, cnr-2000, in http://law.di.unimi.it/webdata/cnr-2000

There you can see the follwoing files available for download:

cnr-2000.graph<br/>
cnr-2000.properties<br/>
cnr-2000-t.graph<br/>
cnr-2000-t.properties

...<br/>
(you can ignore the rest of the files)

The first two files are for the forward (regular) cnr-2000 graph. The other two are for the transpose (inverse) graph. 

What is missing is the "offsets" file. This can be easily created by running:

<pre>
java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000
</pre>
<pre>
java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000-t
</pre>

In our implementations, we assume that the graph deos not have any self-loop. Self-loops can be removed by running:
<pre>
java -Xmx8g -cp "bin":"lib/*" SelfLoopRemover cnr-2000 cnr-2000
</pre>
<pre>
java -Xmx8g -cp "bin":"lib/*" SelfLoopRemover cnr-2000-t cnr-2000-t
</pre>
where the flag "Xmx" specifies the maximum memory allocation pool for a Java virtual machine (JVM). It can be specified in different sizes, such as kilobytes, megabytes, and so on.

Now, we generate weights which are uniformly distibuted. Here, we show this for cnr-2000:
<pre>
java -Xmx8g -cp "bin":"lib/*" GenerateWeightedGraphRandomLong cnr-2000 1 100
</pre>
The above java code produces random weights between range 1 and 100. For each edge, the weight is stored as an integer in the Long format. In our implementations, we access the actual probability of an edge by multplying its corresponding weight by <img src="https://render.githubusercontent.com/render/math?math=10^{-2}">. For instance, for an edge with weight 60, the corresponding probability is obtained by multiplying 60 by 0.01 which is equal to 0.6.

Our implementations work with undirected graphs with symmetrized weights. To change a graph to an undirected one, for each edge we add its inverse. This can be achieved by taking the union of the graph with its transpose. Here, we show how to do this for cnr-2000. Here, we show how to do this for cnr-2000:
<pre>
java -Xmx8g -cp "bin":"lib/*" TransposeWeightedGraphLong cnr-2000 
</pre>
<pre>
java -Xmx8g -cp "bin":"lib/*" SymmetrizeWeightedGraphLong cnr-2000 cnr-2000-t cnr-2000-u
</pre>
The last code creates three files: cnr-2000-u.w.labeloffsets, cnr-2000-u.w.labels, and cnr-2000-u.w.properties. For this dataset, cnr-2000-u.w should be passed as the first input parameter.

Other parametres that are used for ruuning our codes are: threshold, numSamples, and precision.

threshold: It is used to give certainty when outputting the results in a probabilistic graph. 

numSamples: It is used to indicate the number of samples for weakly-global and global nucleus decompositions which are based on Monte-Carlo sampling

precision: It is used to change weights to the actual probabilities. For cnr-2000, precision is equal to 2.

# Compiling

<pre>
mkdir -p bin; javac -cp "bin":"lib/*" -d bin src/it/unimi/dsi/webgraph/labelling/*.java src/*.java
</pre>

# Running
PNucl_DP_Efficient:
<pre>
java -Xmx12g -cp "bin":"lib/*" PNucl_DP_Efficient basename.w precision threshold
</pre>
e.g.
<pre>
java -Xmx12g -cp "bin":"lib/*" PNucl_DP_Efficient newTest.w 16 0.1
</pre>
(Change : to ; if you are on Windows)

The result will be stored in a text file basename+"eta-" + eta + "-bz.txt". The lines of the file are of the form vertex-id:probabilistic core-number.

PNucl_CLT_BinoJava_DP_Efficient_RegularizedBeta:
<pre>
java -Xmx12g -cp "bin":"lib/*" PNucl_CLT_BinoJava_DP_Efficient_RegularizedBeta basename.w precision threshold  
</pre>
e.g.
<pre>
java -Xmx12g -cp "bin":"lib/*" PNucl_CLT_BinoJava_DP_Efficient_RegularizedBeta newTest.w 0.1 1500 2  
</pre>
The result will be stored in a text file basename+"eta-" + eta + "-vc.txt". The lines of the file are of the form vertex-id:probabilistic core-number.

Weakly_Global_onSmallGraph_ForExperiments:
<pre>
java -Xmx12g -cp "bin":"lib/*" Weakly_Global_onSmallGraph_ForExperiments basename.w precision threshold numSamples  
</pre>
e.g.
<pre>
java -Xmx12g -cp "bin":"lib/*" Weakly_Global_onSmallGraph_ForExperiments newTest.w 2 0.1 200 
</pre>
The result will be stored in a text file basename+"eta-" + eta + "-vc.txt". The lines of the file are of the form vertex-id:probabilistic core-number.

Golbal_nuclei_Finding_onSmallGraph_ForExperiments:
<pre>
java -Xmx12g -cp "bin":"lib/*" Golbal_nuclei_Finding_onSmallGraph_ForExperiments basename.w precision threshold numSamples   
</pre>
e.g.
<pre>
java -Xmx12g -cp "bin":"lib/*" Golbal_nuclei_Finding_onSmallGraph_ForExperiments newTest.w 2 0.1 200 
</pre>
The result will be stored in a text file basename+"eta-" + eta + "-vc.txt". The lines of the file are of the form vertex-id:probabilistic core-number.

# Using git
First clone repo.

<pre>
git clone https://github.com/thomouvic/pcore.git
</pre>

This will create a directory "pcore" with the current code of this project. The subdirectories created are "src" and "lib". 

Copy the source files you changed to "pcode/src". 

While being in "pcode", run 
<pre>
git add .
</pre>

Commit changes by running
<pre>
git commit -m "some comment about your changes"
</pre>
If it is the first time, you will be required to run two other commands before. 
<pre>
git config --global user.email "your email"
git config --global user.name "your name"
</pre>

Finally, run
<pre>
git push
</pre>

You will be required to specify username and password. 
If successful, the changes will be in the repository.

