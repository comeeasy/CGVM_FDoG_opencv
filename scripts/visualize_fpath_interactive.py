import sys
import os

import cv2
import numpy as np
import matplotlib.pyplot as plt



def get_paths_from_cli_args():
    # vtf   : value through fpath
    # fpath : fpath
    # target: target sketch img
    
    vtf_path, fpath_path, target_path = sys.argv[1], sys.argv[2], sys.argv[3]
    
    if vtf_path is None:
        print(f"[WARNING] vtf_path is set to default value")
        vtf_path = "/mnt/d/Datasets/simple_fpath_data/fpath_npzs/color_900_fpath_of_infodraw.npz"
    if fpath_path is None:
        print(f"[WARNING] fpath_path is set to default value")
        fpath_path = "/mnt/d/Datasets/simple_fpath_data/fpath_npzs/color_900_fpath.npz"
    if target_path is None:
        print(f"[WARNING] target_path is set to default value")
        target_path = "/home/joono/CGVM_FDoG_opencv/inputs/simple_data/targets/line_900.png"
    
    print(f"[INFO] {vtf_path=}, {fpath_path=}, {target_path=}")
    return vtf_path, fpath_path, target_path

def load_fpath(fpath_path):
    fpath = np.load(fpath_path)["data"]
    print(f"[INFO] Load fpath, {fpath.shape=}")
    return fpath

def load_vtf(vtf_path):
    vtf = np.load(vtf_path)["data"]
    print(f"[INFO] Load vtf, {vtf.shape=}")
    return vtf

def load_target_img(target_path):
    target_img = cv2.imread(target_path) # H x W x 3
    target_img = cv2.cvtColor(target_img, cv2.COLOR_BGR2GRAY) # H x W
    target_img = target_img.transpose(1, 0) # W x H
    print(f"[INFO] Load target img as a grayscale, {target_img.shape=}")
    return target_img

def get_infodraw_from_vtf(vtf):
    # fpath[x, y, 10] is infodraw value
    W, H, _ = vtf.shape
    infodraw = np.zeros((W, H))
    infodraw = vtf[:, :, 10]
    print(f"[INFO] extract infodraw from fpath, {infodraw.shape=}")
    return infodraw
    
# Function to handle mouse events
def handle_mouse_event(event, x, y, flags, param):
    x, y = y, x # We use W x H but opencv uses H x W
    infodraw, target_img, fpath, vtf = param
    
    _, _, threshold_S, _ = fpath.shape
    vtf_box_size = 40
    
    if event == cv2.EVENT_LBUTTONDOWN:
        
        # Draw red line and vtf to left-upper corner through fpath
        for sn in range(threshold_S):
            cur_x, cur_y = fpath[x, y, sn]
            if cur_x > 0 and cur_y > 0:
                px, py = round(cur_x), round(cur_y)
                print(f"{cur_x=}, {cur_y=}, {px=}, {py=}")
                
                # Draw red line
                cv2.rectangle(infodraw, (py, px), (py+1, px+1), color=(0, 0, 255), thickness=1)
                
                # Draw vtf to left-upper corner
                pixel_value = round(vtf[x, y, sn] * 255)
                # fill color
                cv2.rectangle(infodraw, (sn*vtf_box_size, 0), (sn*vtf_box_size+vtf_box_size, 0+vtf_box_size), color=(pixel_value, pixel_value, pixel_value), thickness=-1)
                # border line
                border_color = (0, 0, 255) if sn != 10 else (255, 0, 0)
                cv2.rectangle(infodraw, (sn*vtf_box_size, 0), (sn*vtf_box_size+vtf_box_size, 0+vtf_box_size), color=border_color, thickness=1)
                
        # Draw blue dot to indicate <x, y> for both infodraw and target_img
        cv2.rectangle(infodraw, (y, x), (y+1, x+1), color=(255, 0, 0), thickness=2)
        cv2.rectangle(target_img, (y, x), (y+1, x+1), color=(10, 180, 250), thickness=2)
        
        print(f"Left button clicked at position ({x}, {y})")
    elif event == cv2.EVENT_RBUTTONDOWN:
        print(f"Right button clicked at position ({x}, {y})")    

def display(infodraw, target_img, fpath, vtf):
    
    # cvt [0, 1] -> [0, 255] and BGR colorspace for visualization
    infodraw = np.array(infodraw * 255).astype(np.uint8)
    infodraw = cv2.cvtColor(infodraw, cv2.COLOR_GRAY2BGR)
    target_img = np.array(target_img).astype(np.uint8)
    target_img = cv2.cvtColor(target_img, cv2.COLOR_GRAY2BGR)
    
    # Copy origin to revert to origin img.
    infodraw_origin = infodraw.copy()
    target_img_origin = target_img.copy()
    
    cv2.namedWindow('infodraw')
    cv2.namedWindow('target_img')
    cv2.setMouseCallback('infodraw', handle_mouse_event, (infodraw, target_img, fpath, vtf))
    cv2.setMouseCallback('target_img', handle_mouse_event, (infodraw, target_img, fpath, vtf))
    
    while True:
        
        cv2.setMouseCallback('infodraw', handle_mouse_event, (infodraw, target_img, fpath, vtf))
        cv2.setMouseCallback('target_img', handle_mouse_event, (infodraw, target_img, fpath, vtf))
        cv2.imshow('infodraw', infodraw)
        cv2.imshow('target_img', target_img)
        
        key = cv2.waitKey(1) & 0xFF
        if key == 27:  # ESC key
            print("ESC pressed, exiting...")
            break
        elif key == ord('q'):
            print("'c' key pressed, exiting...")
            break
        elif key == ord('c'):
            infodraw = infodraw_origin.copy()
            target_img = target_img_origin.copy()
    
def main():
    
    vtf_path, fpath_path, target_path = get_paths_from_cli_args()

    fpath = load_fpath(fpath_path)
    vtf = load_vtf(vtf_path) # value through fpath
    target_img = load_target_img(target_path)
    infodraw = get_infodraw_from_vtf(vtf)

    display(infodraw=infodraw, target_img=target_img, fpath=fpath, vtf=vtf)



if __name__ == "__main__":
    main()
    
    