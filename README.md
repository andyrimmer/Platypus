#Platypus

![Build Status](https://github.com/huynguyen250896/Platypus#)

The Platypus variant caller. This version fully supports **Python 3**. The original repository is https://github.com/andyrimmer/Platypus.

---

## Installation
You can install and run Platypus using one of the following two methods:



### Option 1: Using Docker (Recommended)

Using Docker is the easiest way to ensure all dependencies (HTSlib, Cython, Numpy, Pysam) are correctly configured.

1.  **Clone the repository:**
    ```
    git clone [https://github.com/IvantheDugtrio/Platypus.git](https://github.com/IvantheDugtrio/Platypus.git)

    cd Platypus
    ```

2.  **Build the image:**
    ```
    docker build -t platypus:latest .
    ```

3.  **Run Platypus:**

    ```
    docker run -it --rm -v $(pwd):/data platypus:latest python3 Platypus.py --help
    ```

### Option 2: Build from Source

If you prefer to run it natively on your machine, ensure you meet the prerequisites below:

#### System Requirements

* **Python:** 3.10+
* **Compiler:** `gcc` / `g++` 9.0 or greater (build-essential)
* **Libraries:** HTSlib 1.15+ (1.23.1 is recommended)
* **Python packages:** `cython`, `numpy`, `pysam` (Install via `pip install -r requirements.txt`).

#### Build Steps

1.  **Install HTSlib:**

    Download from [HTSlib Web site](http://www.htslib.org/download/).
    ```
    cd htslib-x.x.x

    ./configure && make && sudo make install
    ```
    *Note: If installed in a non-standard location, ensure you export `C_INCLUDE_PATH`, `LIBRARY_PATH`, and `LD_LIBRARY_PATH` accordingly.*

2.  **Build Platypus:**
    ```
    git clone [https://github.com/huynguyen250896/Platypus.git](https://github.com/huynguyen250896/Platypus.git)

    cd Platypus

    pip install -r requirements.txt
    
    make clean

    make PYTHON=python3
    ```
---

## Execution

Once installed, you can run Platypus using:

```

python3 Platypus.py --bamFiles=BAM.bam --refFile=REF.fa --output=variants.vcf

**Note:** For legacy support or older environments (Python 2.7), please refer to the original repository: https://github.com/andyrimmer/Platypus.
