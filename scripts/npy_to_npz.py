import os
import numpy as np

import argparse



def npy_to_npz(input_dir, output_dir):
    # Iterate over all .npy files in the specified directory
    for file in os.listdir(input_dir):
        if file.endswith(".npy"):
            # Load the .npy file
            data = np.load(os.path.join(input_dir, file))
            
            # Define the output path for the .npz file
            output_file = os.path.splitext(file)[0] + ".npz"
            output_path = os.path.join(output_dir, output_file)
            
            # Save as .npz with compression
            np.savez_compressed(output_path, data=data)

            # Optionally, remove the .npy file if no longer needed
            os.remove(os.path.join(input_dir, file))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--input_dir", default="", required=True, help="input directory in which .npy exists")
    parser.add_argument("--output_dir", default="", required=True, help="output directory in which .npz are saved")
    
    args = parser.parse_args()
    
    # Example usage
    npy_to_npz(args.input_dir, args.output_dir)    
