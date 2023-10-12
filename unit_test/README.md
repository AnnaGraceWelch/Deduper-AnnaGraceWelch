This file contains information for the specific cases laid out in my test input SAM file named test_input.sam in this directory. 

My test output file is test_output.sam in this directory. 

There is also a test file containing a list of UMIs called test_umi.txt. 

The first line in my SAM file should be written to output because it cannot be a duplicate if we're retaining first instances. 

The second line is a duplicate exactly the same as the first, so it should NOT be written to output. 

The third line is soft-clipped and the forward strand, so after position adjustment, it should be recognized as a duplicate and NOT written to output. 

The fourth line has a different position, so it SHOULD be written to output.

The fifth line has different strandedness (it's on the minus instead of plus), so it SHOULD be written to output.

The sixth line is soft-clipped and the reverse strand, so after adjusting position, it should be identical to the fifth line. Therefore, it should NOT be written to output. 

The seventh line has a different chromosome than anything we've seen so far, so it SHOULD be written to output.

The eighth line has a different UMI than the previous record, so it SHOULD be written to output. 

Total record lines in input: 8
Total record lines in output: 5