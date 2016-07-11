# Quicksort

This work is using OpenCL description of QuickSort algorithm as a source code for Vivado HLS.

### Introduction to quicksort algorithm:

The quicksort algorithm is developed by Tony Hoare while working in Moscow State University as a vising student, at that time the algorithm was developed during the project of machine translation which was supposed to sort words before looking up to dictionary. His first idea to implement sorting algorithm was insertionsort which was slow and made hime to come up with a new approach that was Quicksort.


Quicksort gained widespread adoption, appearing, for example it became default function library in Unix and C. The basic idea behind Quicksort algorithm is to partition a sequence to be sorted around pivot value, which could be selected in different ways, which all the elements less than pivot value are placed in the left side of the array, all the elements greater than pivot are placed in the right side of the array and all values equal to pivot are placed in the middle of the array.   

![Quicksort] (https://github.com/mediroozmeh/Quicksort/blob/master/Untitled.jpg)


### Quicksort algorithm on GPU:
The nub of the Quicksort algorithm is partitioning a sequence of data to small enough one which can be sorted by single work group, in fact Quicksort algorithm on GPU has two main phases, in the first phase one kernel divide the complete sequence into small enough one which can be sorted by one workgroup, for this purpose the GPU Quicksort function iteratively lunches  **gqsort_kernel** and in the second phase **lqsort_kernel** will complete the sorting job. 





