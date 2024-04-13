import os 
import argparse
import numpy as np 
import scipy.io as sio
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

def generate_RectangularGrid_XZ(Mx, Mz, FoV_azimuth, sensor_width, aspect_ratio):
    FoV_azimuth_radian = FoV_azimuth * (np.pi/180)
    focal_length = (sensor_width/2) / np.tan(FoV_azimuth_radian/2)
    sensor_height = sensor_width / aspect_ratio
    sensor_height_spacing = sensor_height / Mz 
    sensor_width_spacing = sensor_width / Mx

    x = np.arange((sensor_width-sensor_width_spacing)/2, (-sensor_width+sensor_width_spacing)/2 - sensor_width_spacing, -sensor_width_spacing)
    z = np.arange((sensor_height-sensor_height_spacing)/2, (-sensor_height+sensor_height_spacing)/2 - sensor_height_spacing, -sensor_height_spacing)

    X, Z = np.meshgrid(x, z) # Shapes of X and Z: [Mz, Mx]
    xz_grid = np.concatenate((np.expand_dims(X, -1), np.expand_dims(Z, -1)), axis=2)

    FoV_elevation_radian = 2 * np.arctan(sensor_height/(2*focal_length))
    FoV_elevation = (FoV_elevation_radian) / (np.pi/180)

    return xz_grid, focal_length, FoV_elevation

def subtract_background(range_map, background, output_dir):
    # Generate the X-Z rectangular sensing grid
    FoV_azimuth = 100 # degree
    sensor_width = 0.032 # meter
    aspect_ratio = 4/3
    Mx = 160
    Mz = 160
    xz_grid, focal_length, _ = generate_RectangularGrid_XZ(Mx, Mz, FoV_azimuth, sensor_width, aspect_ratio)

    # Zenith angle
    theta_coordinate = (np.pi/2) - np.arctan2(xz_grid[:, :, 1], np.sqrt((xz_grid[:, :, 0])**2 + (focal_length * np.ones((Mz, Mx)))**2 ))

    # Azimuth angle
    phi_coordinate = np.arctan2(focal_length * np.ones((Mz, Mx)), xz_grid[:, :, 0])

    # Estimate the depth maps (DM)
    depth_map = range_map * np.sin(phi_coordinate) * np.sin(theta_coordinate)
    DM_background = background * np.sin(phi_coordinate) * np.sin(theta_coordinate)

    # Background subtraction (BS)
    BS_depth_map = DM_background - depth_map

    # Filter out the undesired reflections
    filtered_BS_depth_map = np.copy(BS_depth_map)
    filtered_BS_depth_map[np.where(filtered_BS_depth_map < 0)] = 0

    # Convert the depth map into binary map
    binary_BS_depth_map = np.copy(filtered_BS_depth_map)
    binary_BS_depth_map[np.where(binary_BS_depth_map > 0)] = 1

    # Prepare DBSCAN
    eps = 2
    min_samples = 1
    DBSCAN_obj = DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean', metric_params=None, algorithm='auto', leaf_size=30, p=None, n_jobs=None)

    idx_z, idx_x = np.where(binary_BS_depth_map==True)
    points = np.concatenate((np.expand_dims(idx_z, axis=1), np.expand_dims(idx_x, axis=1)), axis=1)

    # Clustering
    clustered_depth_map = np.zeros((Mz, Mx))
    clusters = DBSCAN_obj.fit(points)

    label_with_max_number = np.argmax(np.bincount(clusters.labels_))
    points_of_user = []
    for i, label in enumerate(clusters.labels_):
        if label == label_with_max_number:
            clustered_depth_map[points[i][0], points[i][1]] = depth_map[points[i][0], points[i][1]]
            points_of_user.append([points[i][0], points[i][1]])
    points_of_user = np.asarray(points_of_user)
    mean_point_of_user = np.round(np.mean(points_of_user, axis=0)).astype(int)
    print("Zenith angle:", round(theta_coordinate[mean_point_of_user[0], mean_point_of_user[1]] / (np.pi/180), 4))
    print("Azimuth angle:", round(phi_coordinate[mean_point_of_user[0], mean_point_of_user[1]] / (np.pi/180), 4))
    print()

    # Save the depth maps
    output_filename = os.path.join(output_dir, "background_subtraction.mat")
    mdic = {'DM': depth_map,
            'DM_background': DM_background,
            'BS_DM': BS_depth_map,
            'filtered_BS_DM': filtered_BS_depth_map,
            'binary_BS_DM': binary_BS_depth_map,
            'final_DM': clustered_depth_map}
    sio.savemat(output_filename, mdic)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Background subtraction.")
    parser.add_argument(
        "-r", "--range_map", required=True, type=str,
        help="path of the estimated range map")
    parser.add_argument(
        "-b", "--background", required=True, type=str,
        help="path of the background range map")
    args = parser.parse_args()

    # Define output directory
    output_dir = os.path.dirname(args.range_map).replace("range_maps", "BS_results")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load depth maps
    range_map = sio.loadmat(args.range_map)['RDM_est2']
    background = sio.loadmat(args.background)['RDM_est2']

    # Subtract the "single" distance from Feeder to RIS reference
    range_map = range_map - 0.0608 
    background = background - 0.0608

    # Subtract the background from the estimated depth map
    subtract_background(range_map, background, output_dir)

    