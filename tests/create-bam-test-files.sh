#!/bin/bash

# Create the bams directory
mkdir -p data/bams

# Create minimal SAM files for testing BAM merging
cat > data/bams/Sample1-1.sam << 'SAM_EOF'
@HD	VN:1.0	SO:coordinate
@SQ	SN:ref	LN:1000
@RG	ID:Sample1-1	SM:Sample1	LB:lib1	PU:unit1	PL:ION_TORRENT
@PG	ID:samtools	PN:samtools	VN:1.22.1	CL:samtools view -bS
@CO	Test comment line
read1	0	ref	1	60	8M	*	0	0	ATCGATCG	IIIIIIII
read2	0	ref	10	60	8M	*	0	0	GCTAGCTA	IIIIIIII
SAM_EOF

cat > data/bams/Sample1-2.sam << 'SAM_EOF'
@HD	VN:1.0	SO:coordinate
@SQ	SN:ref	LN:1000
@RG	ID:Sample1-2	SM:Sample1	LB:lib1	PU:unit2	PL:ION_TORRENT
@PG	ID:samtools	PN:samtools	VN:1.22.1	CL:samtools view -bS
@CO	Test comment line
read3	0	ref	20	60	8M	*	0	0	TTTTAAAA	IIIIIIII
read4	0	ref	30	60	8M	*	0	0	CCCCGGGG	IIIIIIII
SAM_EOF

cat > data/bams/Sample1-3.sam << 'SAM_EOF'
@HD	VN:1.0	SO:coordinate
@SQ	SN:ref	LN:1000
@RG	ID:Sample1-3	SM:Sample1	LB:lib1	PU:unit3	PL:ION_TORRENT
@PG	ID:samtools	PN:samtools	VN:1.22.1	CL:samtools view -bS
@CO	Test comment line
read5	0	ref	40	60	8M	*	0	0	ACGTACGT	IIIIIIII
SAM_EOF

cat > data/bams/Sample2-1.sam << 'SAM_EOF'
@HD	VN:1.0	SO:coordinate
@SQ	SN:ref	LN:1000
@RG	ID:Sample2-1	SM:Sample2	LB:lib2	PU:unit1	PL:ION_TORRENT
@PG	ID:samtools	PN:samtools	VN:1.22.1	CL:samtools view -bS
@CO	Test comment line
read6	0	ref	50	60	9M	*	0	0	AAAAATTTT	IIIIIIIII
read7	0	ref	60	60	8M	*	0	0	CCCCGGGG	IIIIIIII
SAM_EOF

cat > data/bams/Sample2-2.sam << 'SAM_EOF'
@HD	VN:1.0	SO:coordinate
@SQ	SN:ref	LN:1000
@RG	ID:Sample2-2	SM:Sample2	LB:lib2	PU:unit2	PL:ION_TORRENT
@PG	ID:samtools	PN:samtools	VN:1.22.1	CL:samtools view -bS
@CO	Test comment line
read8	0	ref	70	60	8M	*	0	0	TTTTCCCC	IIIIIIII
SAM_EOF

cat > data/bams/Sample3-1.sam << 'SAM_EOF'
@HD	VN:1.0	SO:coordinate
@SQ	SN:ref	LN:1000
@RG	ID:Sample3-1	SM:Sample3	LB:lib3	PU:unit1	PL:ION_TORRENT
@PG	ID:samtools	PN:samtools	VN:1.22.1	CL:samtools view -bS
@CO	Test comment line
read9	0	ref	80	60	8M	*	0	0	GGGGAAAA	IIIIIIII
read10	0	ref	90	60	8M	*	0	0	TTTTCCCC	IIIIIIII
SAM_EOF

# Convert SAM files to BAM files using samtools with proper sorting
for sam_file in data/bams/*.sam; do
    bam_file="${sam_file%.sam}.bam"
    # Convert to BAM and sort to ensure proper coordinate ordering
    samtools view -bS "$sam_file" | samtools sort -o "$bam_file"
    samtools index "$bam_file"
done

# Clean up SAM files
rm data/bams/*.sam

echo "Test BAM files created successfully!"
echo "Files created:"
ls -la data/bams/
