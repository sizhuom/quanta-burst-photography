# Quanta Burst Photography (SIGGRAPH 2020)

[Sizhuo Ma](https://sizhuoma.netlify.app/), [Shantanu Gupta](https://pages.cs.wisc.edu/~sgupta/), [Arin C. Ulku](https://scholar.google.com/citations?user=ngw0frAAAAAJ), [Claudio Bruschini](https://people.epfl.ch/claudio.bruschini/), [Edoardo Charbon](https://people.epfl.ch/edoardo.charbon), [Mohit Gupta](https://wisionlab.com/people/mohit-gupta/)

[[Project Page](https://wisionlab.com/project/quanta-burst-photography/)] [[Paper](https://wisionlab.com/wp-content/uploads/2020/06/quanta_burst_photography_wision.pdf)] [[Video](https://youtu.be/Mk9XGiND0xE)]

## Overview
This is the official MATLAB implementation of of the paper "Quanta Burst Photography". Single-photon avalanche diodes (SPADs) are an emerging sensor technology capable of detecting individual incident photons, and capturing their time-of-arrival with high timing precision. While these sensors were limited to single-pixel or low-resolution devices in the past, recently, large (up to 1 megapixel) SPAD arrays have been developed. These single-photon cameras (SPCs) are capable of capturing high-speed sequences of binary single-photon images with no read noise. We present quanta burst photography, a computational photography technique that leverages SPCs as passive imaging devices for photography in challenging conditions, including ultra low-light and fast motion. Inspired by recent success of conventional burst photography, we design algorithms that align and merge binary sequences captured by SPCs into intensity images with minimal motion blur and artifacts, high signal-to-noise ratio (SNR), and high dynamic range. We theoretically analyze the SNR and dynamic range of quanta burst photography, and identify the imaging regimes where it provides significant benefits. We demonstrate, via a recently developed SPAD array, that the proposed method is able to generate high-quality images for scenes with challenging lighting, complex geometries, high dynamic range and moving objects. With the ongoing development of SPAD arrays, we envision quanta burst photography will find applications in both consumer and scientific photography.

## Quick Start
1. Run initQBP.m to set up MATLAB paths.
2. Download the test data from [[Google Drive](https://drive.google.com/drive/folders/1SfCtAbMyUWFPkcTjift2LiU8jWqsIeLe?usp=sharing)] [[Dropbox](https://www.dropbox.com/scl/fo/xwp4uu30v5ub2dai6wz1m/h?rlkey=6vh2e3l8zun77k9iawiuwvycp&dl=0)]
2. See ./scripts for how to perform image reconstruction on the test sequences.

## Reference
```
@article{ma_quanta_2020,
    title = “Quanta Burst Photography”,
    author = “Ma, Sizhuo and Gupta, Shantanu and Ulku, Arin C. and Brushini, Claudio and Charbon, Edoardo and Gupta, Mohit”,
    journal = “ACM Transactions on Graphics (TOG)”,
    doi = “10.1145/3386569.3392470”,
    volume = “39”,
    number = “4”,
    year = “2020”,
    month = “7”
    publisher = “ACM”
}
```