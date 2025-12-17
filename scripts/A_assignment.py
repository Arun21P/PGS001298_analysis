from pathlib import Path
import gzip
import shutil
import os
import pandas as pd
import matplotlib.pyplot as plt


# Project root = parent of scripts/
PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"






def unzip_gz_file(gz_path, data_dir):
    """
    Unzips a .gz file and saves the unzipped file ONLY in data directory.
    Returns the unzipped file path.
    """
    gz_path = Path(gz_path)
    data_dir = Path(data_dir)

    txt_filename = gz_path.stem          # removes .gz
    txt_path = data_dir / txt_filename   # FORCE data/

    if not txt_path.exists():
        with gzip.open(gz_path, "rb") as f_in:
            with open(txt_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped file saved to: {txt_path}")
    else:
        print(f"Unzipped file already exists: {txt_path}")

    return txt_path



#A. a. Data Preprocessing.
# """
# Preprocesses the PGS scoring file by removing missing values, filtering to autosomal
# chromosomes (1–22), enforcing correct data types, and retaining GRCh38 harmonized
# positions. The cleaned dataset is saved as a tab-delimited file for downstream analysis.
# """

def preprocess_pgs_data(file_path):
    """
    Loads and preprocesses the full PGS scoring file.
    """
    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#"
    )

    df.isna().sum()

    # output
    # rsID                    118
    # chr_name                  0
    # chr_position              0
    # effect_allele             0
    # other_allele              0
    # effect_weight             0
    # hm_source                 0
    # hm_rsID                 120
    # hm_chr                    3
    # hm_pos                    3
    # hm_inferOtherAllele    9227
    # dtype: int64



    df = df.dropna(axis=1, how="all")
    df = df.dropna()
    df.isna().sum()


    # output
    # rsID             0
    # chr_name         0
    # chr_position     0
    # effect_allele    0
    # other_allele     0
    # effect_weight    0
    # hm_source        0
    # hm_rsID          0
    # hm_chr           0
    # hm_pos           0
    # dtype: int64


    # Keep all  chromosomes (autosomes + sex)
    # chr_map = {"X": 23, "Y": 24}

    # df["hm_chr"] = df["hm_chr"].replace(chr_map)
    # df["hm_chr"] = df["hm_chr"].astype(int)


    # Keep only numeric chromosomes (autosomes 1-22)
    df = df[df["hm_chr"].isin([str(i) for i in range(1, 23)])]

    # Correct data types
    df["hm_chr"] = df["hm_chr"].astype(int)
    df["hm_pos"] = df["hm_pos"].astype(int)
    df["effect_weight"] = df["effect_weight"].astype(float)

    # rsID as string
    df["rsID"] = df["rsID"].astype(str)

    grch37_cols = ["chr_name", "chr_position"]
    df = df.drop(columns=grch37_cols)

    print(df.head())
    # Define new file name
    output_file = OUTPUT_DIR/"PGS001298_hmPOS_GRCh38_cleaned.txt"

    # Save as tab-delimited text
    df.to_csv(output_file, sep="\t", index=False)

    print(f"Cleaned DataFrame saved to {output_file}")

    
    return df




# B. b.Exploratory Analysis.
# """
# Performs basic exploratory analysis on the cleaned PGS dataset by reporting dataset size,
# chromosome-wise variant counts, and summary statistics, and visualizing the distribution
# of effect weights across all autosomal chromosomes.
# """

def exploratory_analysis(clean_df):

    # Number of variants and columns
    print("Number of variants:", clean_df.shape[0])
    print("Number of columns:", clean_df.shape[1])

    # Summary stats for numeric columns
    clean_df[["hm_chr", "hm_pos", "effect_weight"]].describe()

    # Unique chromosomes
    print("Unique chromosomes in dataset:", sorted(clean_df["hm_chr"].unique()))

    # Count of variants per chromosome
    print("\nVariants per chromosome:")
    print(clean_df["hm_chr"].value_counts().sort_index())



    plt.figure(figsize=(6,4))
    plt.hist(clean_df["effect_weight"], bins=50)
    plt.xlabel("Effect Weight")
    plt.ylabel("Frequency")
    plt.title("Distribution of Effect Weight (All Chromosomes)")
    plt.show()
 # create folder if missing
    save_path = OUTPUT_DIR / f"Distribution_of_Effect_Weight.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Distribution_of_Effect_Weight Plot saved to {save_path}")





# C. c. Generate a histogram to depict the distribution of the column “effect_weight” for all the locations of chromosome 21. Only the columns “hm_chr” and “hm_pos” need to be considered for chromosome locations.
# """
# Generates and saves a histogram of effect weights for a specified chromosome,
# summarizing the distribution of variant effect sizes and reporting basic
# descriptive statistics for that chromosome.
# """
def plot_chr_effect_weight_histogram(clean_df, chr_num=None):
    """
    Plots the distribution of effect_weight for a given chromosome
    and saves the plot in the specified directory.
    
    Parameters:
        clean_df (pd.DataFrame): Cleaned PGS DataFrame.
        chr_num (int): Chromosome number to filter).
    """
    
    # Filter the specified chromosome
    chr_df = clean_df[clean_df["hm_chr"] == chr_num][["hm_chr", "hm_pos", "effect_weight"]]

    if chr_df.empty:
        print(f"No variants found for chromosome {chr_num}.")
        return

    # Plot histogram
    plt.figure(figsize=(6, 4))
    plt.hist(chr_df["effect_weight"], bins=50)
    plt.xlabel("Effect Weight")
    plt.ylabel("Frequency")
    plt.title(f"Effect Weight Distribution on Chromosome {chr_num}")
    plt.show()



 # create folder if missing
    save_path = OUTPUT_DIR / f"chr{chr_num}_effect_weight_hist.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Effect Weight Distribution on Chromosome {chr_num} Plot saved to {save_path}")

    # Print variant stats
    print(f"Number of variants on chromosome {chr_num}: {chr_df.shape[0]}")
    print(chr_df["effect_weight"].describe())


if __name__ == "__main__":
    
    # Step 1: Unzip

    gz_file = DATA_DIR / "PGS001298_hmPOS_GRCh38.txt.gz"
    txt_file = unzip_gz_file(gz_file, DATA_DIR)

    # OR 
    # clean_df = pd.read_csv(gz_file, sep="\t", comment="#", compression='gzip')

    
    # # Step 2: Preprocessing (full file)
    clean_df = preprocess_pgs_data(txt_file)
    
    # # Step 3: Exploratory analysis (full file)
    exploratory_analysis(clean_df)
    
    # # Step 4: Chromosome 21 histogram
    plot_chr_effect_weight_histogram(clean_df, chr_num=21)

    print("PGS001298 Analysis Completed")
    print("Conclusion:The effect weight distribution across all chromosomes is approximately symmetric and centered around zero, indicating that most variants have small additive effects, as expected for a polygenic trait. Chromosome 21 shows a similar pattern but with fewer variants and a slightly narrower spread, suggesting no chromosome-specific bias and a consistent contribution of chr21 variants to the overall polygenic score.")




