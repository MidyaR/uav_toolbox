"""
                     Metashape residual error exporter                 
This script exports the residual error for each tie point in a Metashape project.
The exported file contains the tie point ID, the X, Y, Z coordinates in the project's CRS,
the projection coordinates in pixel units, and the root-mean-square error (RMSE) of the tie point.
The RMSE is calculated as the square root of the average squared reprojection error of the tie point
in all the images where it is visible.
Note:
- This script assumes that the tie points in the chunk are valid and have been optimized.
- This script does not export the error of control points or checkerboard targets.
Information:
Author:         Midya Rostami
Partners: Marjan Ahangarha, Masood Varshosaz
Author_Email:   midyalab@gmail.com 
Author_Website: https://github.com/MidyaR/uav_toolbox/
Author_Address: ​Close Range Photogrammetry & Robotics Lab, Department of Photogrammetry and Remote Sensing, 
                Faculty of Geodesy and Geomatics Engineering, K. N. Toosi University of Technology, Tehran, Iran.19967-15433

********************** License: All Rigths Reserved ©2023  ********************
"""
import Metashape
import math, time, statistics
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
# import pandas as pd

# Checking compatibility
compatible_major_version = "1.8"
found_major_version = ".".join(Metashape.app.version.split('.')[:2])
if found_major_version != compatible_major_version:
    raise Exception("Incompatible Metashape version: {} != {}".format(found_major_version, compatible_major_version))

def export_point_cloud():
    doc = Metashape.app.document
    chunk = doc.chunk
    M = chunk.transform.matrix
    crs = chunk.crs
    point_cloud = chunk.point_cloud
    projections = point_cloud.projections
    points = point_cloud.points
    npoints = len(points)
    tracks = point_cloud.tracks

    path = Metashape.app.getSaveFileName("Specify export file name:", filter=" *.txt")
    if not path:
        print("Incorrect path, script aborted.")
        return None	
    file = open(path, "wt")
    print("Script started")

    t1 = time.time()
    points_errors = {}

    for photo in chunk.cameras:
        if not photo.transform:
            continue

        T = photo.transform.inv()
        calib = photo.sensor.calibration

        point_index = 0
        for proj in projections[photo]:
            track_id = proj.track_id
            while point_index < npoints and points[point_index].track_id < track_id:
                point_index += 1
            if point_index < npoints and points[point_index].track_id == track_id:
                if not points[point_index].valid:
                    continue

                coord = T * points[point_index].coord
                coord.size = 3
                dist = calib.error(coord, proj.coord).norm() ** 2
                v = M * points[point_index].coord
                v.size = 3

                if point_index in points_errors.keys():
                    point_index = int(point_index)
                    points_errors[point_index].x += dist
                    points_errors[point_index].y += 1
                else:
                    points_errors[point_index] = Metashape.Vector([dist, 1])

    total_error = 0
    valid_points = 0
    errors = []
    X_values = []
    Y_values = []

    for point_index in range(npoints):
        if not points[point_index].valid:
            continue

        if chunk.crs:
            w = M * points[point_index].coord
            w.size = 3
            X, Y, Z = chunk.crs.project(w)
        else:
            X, Y, Z, w = M * points[point_index].coord

        X_values.append(X)
        Y_values.append(Y)

        error = math.sqrt(points_errors[point_index].x / points_errors[point_index].y)

        errors.append(error)
        total_error += error#
        valid_points += 1#

        file.write("{:6d}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(point_index, X, Y, Z, error))


    mean_error = total_error / valid_points
    std_dev = statistics.stdev(errors)
    skewness = stats.skew(errors)
    kurtosis = stats.kurtosis(errors)
    print("Mean Error: ", mean_error, "px")
    print("Standard Deviation: ", std_dev)
    print("Skewness: ", skewness)
    print("Kurtosis: ", kurtosis)
    

    #plots 
    fontsizes= 16
    fig,ax = plt.subplots(1)
    plt.scatter(X_values, Y_values, s=1, c=errors, cmap='jet')
    plt.colorbar(label='Residual Classes')
    plt.xlabel('X', fontsize=fontsizes)
    plt.ylabel('Y', fontsize=fontsizes)
    plt.title('Residual Error Map (px)', fontsize=fontsizes, fontweight='bold')
    plt.savefig('residual_map.png',dpi = 800)
    plt.show()

    fig,ax = plt.subplots(1)
    plt.hist(errors, bins='auto', edgecolor='black')
    plt.xlabel("Values (px)", fontsize=fontsizes)
    plt.ylabel("Frequency", fontsize=fontsizes)
    plt.title("Distribution of Residuals", fontsize=fontsizes, fontweight='bold')
    plt.grid(False)
    plt.savefig('residual_histogram.png',dpi = 800)
    plt.show()


    t2 = time.time() 
    file.flush()
    file.close()

    # Save the parametrs  
    path2 = Metashape.app.getSaveFileName("Specify export file name for statisics:", filter=" *.txt")
    if not path:
        print("Incorrect path, script aborted.")
        return None	
    file2 = open(path2, "wt")
    str_num0 = 'Residual error mean: ' + format(std_dev, ".4f")+' px'
    str_num1 = '\nStandard Deviation: ' + format(mean_error, ".4f")
    str_num2 = '\nSkewness: ' + format(skewness, ".4f")
    str_num3 = '\nkurtosis: ' + format(kurtosis, ".4f")
    file2.write(str_num0)
    file2.write(str_num1)
    file2.write(str_num2)
    file2.write(str_num3)
    file2.flush()
    file2.close()

    print("Script finished in " + "{:.2f}".format(t2-t1) + " seconds")
    return True


label = "Export Rsidual Stat"
Metashape.app.addMenuItem(label, export_point_cloud)
print("To execute this script press {}".format(label))
