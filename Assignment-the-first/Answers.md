# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | R1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | I1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | I2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | R2 | 101 | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
![image](https://user-images.githubusercontent.com/106117735/181919670-c4fa3fc8-ccd4-4dcc-b524-e6e0208511a7.png)
![image](https://user-images.githubusercontent.com/106117735/181919681-d78ea8fa-d492-42c2-812e-9994e3982fb6.png)
![image](https://user-images.githubusercontent.com/106117735/181919689-d0b1933c-75cd-461c-957f-f534c5c39083.png)
![image](https://user-images.githubusercontent.com/106117735/181919701-a998caba-f884-4630-885f-bbb8ae79d96e.png)

3. Good quality score cutoff? 
Index files, 32. 
Bio files, 35. Determined by distribution, mostly uniform above the threshold. 
4. Command to count indexes with N character in one line. 

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz| awk 'NR % 2==0'| grep -c "N"``` 

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz| awk 'NR % 2==0'| grep -c "N"```
Returned: 3976613 + 3328051 indexes with "N" = 7304664


## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
