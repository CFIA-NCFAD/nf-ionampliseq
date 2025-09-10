# Create the fastq directory
mkdir -p data/fastq

# Create minimal FASTQ files for Sample1 (3 files)
cat > data/fastq/Sample1-1.fastq << 'EOF'
@Sample1_1
ATCGATCG
+
IIIIIIII
@Sample1_2
GCTAGCTA
+
IIIIIIII
EOF

cat > data/fastq/Sample1-2.fastq << 'EOF'
@Sample1_3
TTTTAAAA
+
IIIIIIII
@Sample1_4
CCCCGGGG
+
IIIIIIII
EOF

cat > data/fastq/Sample1-3.fastq << 'EOF'
@Sample1_5
ACGTACGT
+
IIIIIIII
EOF

# Create minimal FASTQ files for Sample2 (2 files)
cat > data/fastq/Sample2-1.fastq << 'EOF'
@Sample2_1
AAAAATTTT
+
IIIIIIIII
@Sample2_2
CCCCGGGG
+
IIIIIIII
EOF

cat > data/fastq/Sample2-2.fastq << 'EOF'
@Sample2_3
TTTTCCCC
+
IIIIIIII
EOF

# Create minimal FASTQ file for Sample3 (1 file)
cat > data/fastq/Sample3-1.fastq << 'EOF'
@Sample3_1
GGGGAAAA
+
IIIIIIII
@Sample3_2
TTTTCCCC
+
IIIIIIII
EOF

# Compress all files in parallel
parallel gzip ::: data/fastq/*.fastq
