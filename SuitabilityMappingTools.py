"""
The SuitabilityMappingTools.py script is a Python-based geospatial analysis tool designed to process, analyze,
and conduct sensitivity assessments on raster datasets (GeoTIFF files).

The script is specifically built to facilitate weighted overlay analysis for environmental and hydrological
suitability mapping, using raster layers representing different thematic variables
 (e.g., soil properties, permeability, vegetation indices).

"""




import os
import glob
import rasterio
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rasterio.warp import reproject
from rasterio.enums import Resampling
import copy

# Functions---------------------------------
def plot_layers(data_dir):
    # Get list of all raster files in the data directory
    raster_files = [f for f in os.listdir(data_dir) if f.endswith('.tif')]

    # Set up the plot
    num_files = len(raster_files)
    if num_files == 0:
        print("No raster files found in the directory.")
        return

        # Set up the plot
    if num_files > 3:
        cols = int(np.ceil(np.sqrt(num_files)))  # Number of columns
        rows = int(np.ceil(num_files / cols))  # Number of rows
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))
        axes = axes.flatten()  # Flatten for easier indexing
    elif num_files > 1:
        fig, axes = plt.subplots(1, num_files, figsize=(num_files * 5, 5))
    else:
        fig, axes = plt.subplots(1, 1, figsize=(5, 5))
        axes = [axes]  # Wrap single axis in a list for consistent handling

    for ax, raster_file in zip(axes, raster_files):
        with rasterio.open(os.path.join(data_dir, raster_file)) as src:
            data = src.read(1)
            # Mask no-data values
            data = data.astype(float)
            data[data == src.nodata] = float('nan')

            # Determine the colormap based on the filename
            if raster_file.endswith(('_mean.tif', '_cv.tif', '_range.tif', '_std.tif')):
                cmap = 'RdBu_r'  # Red to blue with red reflecting high values and blue reflecting low values
            else:
                cmap = 'BrBG'  # Default colormap

            # Plot the data with scaling based on min and max values
            img = ax.imshow(data, cmap=cmap, vmin=np.nanmin(data), vmax=np.nanmax(data))
            ax.set_title(os.path.basename(raster_file))
            fig.colorbar(img, ax=ax, orientation='vertical')

    plt.tight_layout()
    plt.show()
#-------------------------------
def get_layer_info(filepath):
    """
    This function retrieves information about a TIFF file at a specified path.
    for example: filepath = r'G:\\Data\\Thinning_Recharge_Final_Data\\SensitivityAnalysis\\sBA.tif'
    """

    with rasterio.open(filepath) as src:
        data = src.read(1).astype(np.float32)  # Convert to float32 to avoid data type issues

        # Mask NoData values
        nodata_value = src.nodata
        if nodata_value is not None:
            data = np.where(data == nodata_value, np.nan, data)

        # Compute stats while ignoring NoData values
        min_value = np.nanmin(data)
        max_value = np.nanmax(data)
        mean_value = np.nanmean(data)

        info = {
            'File': os.path.basename(filepath),
            'Filepath': filepath,
            'Projection': src.crs,
            'Min Value': min_value,
            'Max Value': max_value,
            'Mean Value': mean_value,
            'Extent': src.bounds,
            'Width': src.width,
            'Height': src.height,
            'Dimensions': src.shape,
            'Data Type': src.dtypes[0]
        }
        #print(info)
        return info
#---------------------------------
def get_dir_info(tiff_dir):
    """
    This function retrieves information about the TIFF files in the specified directory.
    """
    tiff_files = glob.glob(os.path.join(tiff_dir, '*.tif'))
    tiff_info = []
    for tiff_file in tiff_files:
        info = get_layer_info(tiff_file)
        print(info)
#-------------------------------
def resample_layers(tiff_dir, resampling_method=Resampling.nearest):
    """
    This function resamples the TIFF files in the specified directory to match the dimensions of the first file.
    """

    tiff_files = glob.glob(os.path.join(tiff_dir, '*.tif'))
    if not tiff_files:
        print("No TIFF files found in the directory.")
        return

    # Create resampled directory
    resampled_dir = os.path.join(tiff_dir, 'resampled')
    if not os.path.exists(resampled_dir):
        os.makedirs(resampled_dir)
    else:
        overwrite = input("Resampled directory already exists. Do you want to overwrite the files? (Y/N): ").strip().upper()
        if overwrite != 'Y':
            print("Resampling aborted.")
            return

    # Use the first file as the reference
    with rasterio.open(tiff_files[0]) as ref_src:
        ref_meta = ref_src.meta.copy()
        ref_transform = ref_src.transform
        ref_crs = ref_src.crs
        ref_width = ref_src.width
        ref_height = ref_src.height

    for tiff_file in tiff_files:
        #with rasterio.open(tiff_file) as src:
        #    data = src.read(1)
        with rasterio.open(tiff_file) as src:
            data = src.read(1).astype(np.float32)
            # Handle NoData values
            nodata_value = src.nodata
            if nodata_value is not None:
                data = np.where(data == nodata_value, np.nan, data)

            transform, width, height = rasterio.warp.calculate_default_transform(
                src.crs, ref_crs, ref_width, ref_height, *src.bounds)

            kwargs = src.meta.copy()
            kwargs.update({
                'crs': ref_crs,
                'transform': ref_transform,
                'width': ref_width,
                'height': ref_height,
                'dtype':rasterio.float32,
                'nodata': np.nan
            })

            output_file = os.path.join(resampled_dir, f"{os.path.splitext(os.path.basename(tiff_file))[0]}_rsmp.tif")
            with rasterio.open(output_file, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=ref_transform,
                        dst_crs=ref_crs,
                        resampling=resampling_method
                    )
    print(f"Resampled layers saved in {resampled_dir}")
#---------------------------------
def select_tif_files(directory):
    """
    Lists .tif files in a directory, prompts the user to select files by index,
    retrieves their metadata using get_layer_info, and adds a weight key.

    Args:
        directory (str): Path to the directory containing .tif files.

    Returns:
        list of dict: List of dictionaries containing TIFF metadata and weight.
    """

    # Get all .tif files in the directory
    tif_files = glob.glob(os.path.join(directory, "*.tif"))

    if not tif_files:
        print("No .tif files found in the directory.")
        return []

    # Extract basenames
    basenames = [os.path.basename(f)[:-4] for f in tif_files]  # Remove ".tif"

    # Display available files with index
    print("\nAvailable .tif files:")
    for i, name in enumerate(basenames):
        print(f"{i}: {name}")

    print('Enter the indices of the files you want to select (comma-separated)')

    # Prompt user for selection
    selected_indices = input("\nEnter here:  ")

    try:
        selected_indices = list(map(int, selected_indices.split(",")))
    except ValueError:
        print("Invalid input. Please enter numbers separated by commas.")
        return []

    # Filter selected files
    selected_files = [tif_files[i] for i in selected_indices if i in range(len(tif_files))]

    if not selected_files:
        print("No valid selections made.")
        return []

    # Compute equal weight for all selected layers
    weight = 1 / len(selected_files)

    # Get metadata and add weight to each selected file
    files_equal_weighting = [get_layer_info(f) for f in selected_files]

    for metadata in files_equal_weighting:
        metadata['Weight'] = weight

    return files_equal_weighting
#-------------------------------------------
def define_weights(files_list_dic, user_input=True):
    """
    Defines weights for a list of TIFF files.

    Args:
        files_list_dic (list): List of dictionaries containing file metadata.
        user_input (bool): If True, prompts user for input. If False, assigns equal weights.

    Returns:
        list: Updated list with assigned weights, or None if invalid input.
    """

    # ✅ Auto-assign equal weights if `user_input=False`
    if not user_input:
        num_files = len(files_list_dic)
        weight = 1 / num_files  # Assign equal weight

        for file_info in files_list_dic:
            file_info["Weight"] = weight

        return files_list_dic  # ✅ Return updated list

    # ✅ If `user_input=True`, ask the user to modify weights
    print("\nCurrent file weights:")
    for file_info in files_list_dic:
        print(f"{file_info['File']}: {file_info['Weight']}")

    modify_weights = input("\nDo you want to modify weights? (Y/N): ").strip().upper()

    if modify_weights == "Y":
        print("\nEnter new weights (comma-separated). Make sure they sum to 1.")
        new_weights = input("\nEnter new weights: ").split(",")

        try:
            new_weights = list(map(float, new_weights))
        except ValueError:
            print("❌ ERROR: Invalid weight input. Please enter numbers.")
            return None

        if abs(sum(new_weights) - 1.0) > 1e-6:
            print("❌ ERROR: Weights must sum to 1. Restarting...")
            return None

        # ✅ Assign new user-defined weights
        for i, file_info in enumerate(files_list_dic):
            file_info["Weight"] = new_weights[i]

    return files_list_dic  # ✅ Return updated list
#------------------------------------------
def check_weights(files_list_dic):
    print("\nCurrent file weights:")
    for file_info in files_list_dic:
        print(f"{file_info['File']}: {file_info['Weight']}")

    print("If weights are incorrect, run define_wights() function to change")
#---------------------------------------
def plot_histogram(filepath, title):
    """
    Plots a histogram of raster values, excluding no-data values.

    Args:
        filepath (str): Path to the raster file.
        title (str): Descriptive label for the plot.
    """
    with rasterio.open(filepath) as src:
        data = src.read(1).astype(np.float32)
        nodata = src.nodata
        
        # Handle NoData values properly
        data = np.where(data == -9999, np.nan, data)  # Treat -9999 as NoData
        if nodata is not None:
            data = np.where(data == nodata, np.nan, data)

        # Remove small values (values < 0.9999)
        data = np.where(data < 0.9999, np.nan, data)

        # Remove NaN values BEFORE computing histogram
        data = data[~np.isnan(data)]

        # Ensure there are valid data points left
        if data.size == 0:
            print(f"❌ No valid data to plot for {os.path.basename(filepath)}.")
            return

        # Define bins for plotting
        bins = np.arange(0, 11, 1)

        # Plot histogram
        fig, ax = plt.subplots()
        counts, bins, patches = ax.hist(data, bins=bins, edgecolor='black')

        # Add vertical line for mean value
        nanmean = np.nanmean(data)  # Use nanmean to handle NaN values
        ax.axvline(nanmean, color='red', linestyle='--', label=f'Mean: {nanmean:.2f}')
        ax.legend()

        # Format axes and labels
        ax.set_title(f"Histogram of {os.path.basename(filepath)}: {title}")
        ax.set_xlabel("Value")
        ax.set_ylabel("Frequency (thousands of pixels)")
        ax.set_xticks(bins)
        ax.set_yticklabels([f'{int(label/1000)}k' for label in ax.get_yticks()])

        # Count of valid pixels (excluding NoData)
        pixel_count = data.size

        # Calculate percentages of pixels in each bin
        percentages = {int(bins[i]): (counts[i] / pixel_count) * 100 for i in range(len(counts))}

        # Save the figure
        output_file = os.path.splitext(filepath)[0] + '.png'
        plt.savefig(output_file)
        plt.close()

        # Save the markdown table
        table_file = os.path.splitext(filepath)[0] + '_histogram.txt'
        with open(table_file, 'w') as f:
            f.write(f"# Histogram Data for {os.path.basename(filepath)}\n")
            f.write(f"Total count of pixels (excluding no-data): {pixel_count}\n\n")
            f.write("| Bin | Count | Percentage of Total Pixels | Min Value | Max Value | Area (ha) |\n")
            f.write("|-----|-------|----------------------------|-----------|-----------|-----------|\n")
            for i in range(len(counts)):
                bin_min = bins[i]
                bin_max = bins[i + 1]
                pixels = int(counts[i])
                sqmet = pixels * 900
                hectares = sqmet / 10000
                f.write(f"| {int(bins[i])} | {pixels} | {percentages[int(bins[i])]:.2f}% | {bin_min:.2f} | {bin_max:.2f} | {hectares:.2f} |\n")

        print(f"✅ Histogram saved to {output_file}")
        print(f"✅ Histogram data saved to {table_file}")

        # ✅ Debug: Print out the histogram stats
        print(f"Count of pixels (excluding no-data pixels): {pixel_count}")
        for i in range(len(counts)):
            print(f"Percentage of pixels in bin {int(bins[i])}: count = {int(counts[i])} : {percentages[int(bins[i])]:.2f}%")


#-----------------------------------------
def create_weighted_index(files_list_dic, output_file_name, output_file_directory=False, save_hist=True):
    """
    Creates a weighted index from a list of TIFF files and outputs a weighted raster.
    Ensures NoData handling is consistent with `get_layer_info()`.
    """

    # Set default output directory
    file_info_dict = files_list_dic[0]
    fp = file_info_dict['Filepath']
    dirr = os.path.dirname(fp)

    if not output_file_directory:
        output_file_directory = os.path.join(dirr, "weighted_layers")
        print(f'output_file_directory: {output_file_directory}')

    os.makedirs(output_file_directory, exist_ok=True)

    # Ensure the filename does not have an extra .tif
    if not output_file_name.lower().endswith(".tif"):
        output_file_name += ".tif"

    output_raster_path = os.path.join(output_file_directory, output_file_name)

    if os.path.exists(output_raster_path):
        print(f"Deleting existing file: {output_raster_path}")
        os.remove(output_raster_path)  # Ensure clean overwrite

    print("\nCurrent file weights:")
    wts = []
    for file_info in files_list_dic:
        print(f"{file_info['File']}: {file_info['Weight']}")
        wts.append(file_info['Weight'])

    print("\nGenerating weighted raster...")
    weighted_sum = None

    for file_info in files_list_dic:
        weight = file_info['Weight']
        with rasterio.open(file_info['Filepath']) as src:
            data = src.read(1).astype(np.float32)
            nodata_value = src.nodata
            if nodata_value is not None:
                data = np.where(data == nodata_value, np.nan, data)

        if weighted_sum is None:
            weighted_sum = np.zeros_like(data, dtype=np.float32)

        weighted_sum += np.nan_to_num(data) * weight

    weighted_sum[weighted_sum < 0.999] = np.nan

    # ✅ Ensure NoData handling before writing
    if np.isnan(weighted_sum).all():
        print("❌ ERROR: weighted_sum contains only NaN values. The file will not be valid.")
        return None  # Exit early

    weighted_sum = np.nan_to_num(weighted_sum, nan=-9999)  # Convert NaNs to -9999

    # ✅ Debug: Check for NaN, Inf, or all-zero values before writing
    if np.isnan(weighted_sum).any():
        print("❌ ERROR: weighted_sum contains NaN values before saving.")
        weighted_sum = np.nan_to_num(weighted_sum)

    if np.isinf(weighted_sum).any():
        print("❌ ERROR: weighted_sum contains infinite values.")
        weighted_sum = np.nan_to_num(weighted_sum)

    if np.all(weighted_sum == 0):
        print("❌ ERROR: weighted_sum contains only zeros. The output will be invalid.")
        return None

    # Get metadata and update it
    with rasterio.open(files_list_dic[0]['Filepath']) as src:
        meta = src.meta.copy()
    meta.update(dtype=rasterio.float32, nodata=-9999)  # ✅ Fix: Use -9999 instead of NaN for NoData

    # ✅ Debug: Print metadata before writing
    #print("Raster Metadata Before Writing:")
    #print(meta)

    # ✅ Save the weighted raster
    try:
        with rasterio.open(output_raster_path, 'w', **meta) as dst:
            dst.write(weighted_sum, 1)

        print(f"✅ Weighted raster saved successfully: {output_raster_path}")
    except rasterio.errors.RasterioIOError:
        print(f"WARNING: Could not write to {output_raster_path}. Close open instances and try again.")
        return None

    # ✅ Debug: Ensure the file was written correctly
    if os.path.exists(output_raster_path):
        file_size = os.path.getsize(output_raster_path)
        print(f"✅ File {output_raster_path} saved successfully. Size: {file_size} bytes")
    else:
        print(f"❌ ERROR: {output_raster_path} was not created.")
        return None

    # ✅ Save metadata to a separate text file
    output_metadata_path = output_raster_path.replace('.tif', '_metadata.txt')
    with open(output_metadata_path, 'w') as f:
        f.write(f"Date: {pd.Timestamp.now()}\nFiles combined:\n")
        for file_info in files_list_dic:
            f.write(f"{file_info['File']} with weight {file_info['Weight']}\n")

    print("\nMetadata saved successfully!")

    # ✅ Ensure the file is fully closed before reopening
    del dst

    # ✅ Generate histogram if requested
    if save_hist:
        plot_histogram(output_raster_path,wts)

    return output_raster_path  # Return the created file path
#--------------------------------------
def sensitivity_analysis_weights(files_list_dic, run_name):
    """
    Conducts sensitivity analysis by shuffling weights and computing weighted suitability maps.
    """

    sensitivity_dir = os.path.join(os.path.dirname(files_list_dic[0]['Filepath']), "Sensitivity_Analysis", run_name)
    os.makedirs(sensitivity_dir, exist_ok=True)

    file_names = [file_info['File'] for file_info in files_list_dic]
    base_weights = [file_info['Weight'] for file_info in files_list_dic]
    weight_permutations = list(itertools.permutations(base_weights))

    results = []
    output_files = []

    for i, weight_set in enumerate(weight_permutations):
        for j, file_info in enumerate(files_list_dic):
            file_info["Weight"] = weight_set[j]

        output_filename = f"{run_name}_WSA_{i}"

        # 🛠 **Get output path from `create_weighted_index()`**
        output_path = create_weighted_index(files_list_dic, output_filename, sensitivity_dir)

        # **Ensure the file was created before proceeding**
        if output_path is None or not os.path.exists(output_path):
            print(f"❌ Error: Expected output file {output_filename} was not created.")
            continue  # Skip this iteration if the file wasn't created

        # Save path to use for creating summary stats
        output_files.append(output_path)

        # Analyze the newly created raster
        with rasterio.open(output_path) as src:
            data = src.read(1)
            nodata = src.nodata

        if nodata is not None:
            data = data[data != nodata]  # Remove NoData values
            if data.size == 0:
                print(f"❌ Warning: No valid data in {output_filename}. Skipping.")
                continue

        percent_above_5 = np.mean(data > 5) * 100
        percent_above_6 = np.mean(data > 6) * 100
        percent_above_7 = np.mean(data > 7) * 100
        d = data
        d[d == -9999] = np.nan
        nanmean_value = np.nanmean(d)

        results.append([output_filename] + list(weight_set) + [percent_above_5, percent_above_6, percent_above_7, nanmean_value])


    results_df = pd.DataFrame(results, columns=["File", *file_names, "% > 5", "% > 6", "% > 7", "Mean"])
    results_csv_path = os.path.join(sensitivity_dir, f"{run_name}_WSA_combined_results.csv")
    results_df.to_csv(results_csv_path, index=False)

    print(f"✅ Sensitivity analysis completed. Results saved at: {results_csv_path}")
    return output_files, sensitivity_dir, run_name
#------------------------------
def sensitivity_analysis_removal(files_list_dic, run_name):
    """
    Conducts sensitivity analysis by removing each input variable one at a time and reweighting.
    """
    remov = '_'
    # ✅ Ensure valid file list
    if not files_list_dic or len(files_list_dic) == 0:
        print("❌ ERROR: No valid files provided for analysis. Exiting.")
        return

    # ✅ Ensure Sensitivity_Analysis directory exists
    sensitivity_dir = os.path.join(os.path.dirname(files_list_dic[0]['Filepath']), "Sensitivity_Analysis", run_name)
    os.makedirs(sensitivity_dir, exist_ok=True)

    # ✅ Control scenario with all variables equally weighted
    num_files = len(files_list_dic)

    files_list_copy = copy.deepcopy(files_list_dic)

    for j, file_info in enumerate(files_list_copy):
        file_info['Weight'] = 1 / num_files


    control_output_filename = f"{run_name}_RM_control.tif"
    print(f"Generating control scenario: {control_output_filename}")

    control_output_path = create_weighted_index(files_list_dic, control_output_filename, sensitivity_dir)

    if control_output_path is None or not os.path.exists(control_output_path):
        print(f"❌ ERROR: Expected control file {control_output_filename} was not created. Skipping analysis.")
        return

    plot_histogram(filepath=control_output_path, weights="Control")

    # ✅ Read and analyze control raster
    with rasterio.open(control_output_path) as src:
        data = src.read(1)
        data = data[~np.isnan(data)]
        if data.size == 0:
            print(f"❌ Warning: No valid data in {control_output_path}. Skipping.")
            return

        control_percent_above_5 = np.mean(data > 5) * 100
        control_percent_above_6 = np.mean(data > 6) * 100
        control_percent_above_7 = np.mean(data > 7) * 100
        d = data
        d[d == -9999] = np.nan
        control_nanmean_value = np.nanmean(d)

    # ✅ Store control results
    results = []  # Initialize empty list for storing all results
    results.append(
        ["Control", control_percent_above_5, control_percent_above_6, control_percent_above_7, control_nanmean_value])

    # ✅ Track created rasters
    output_files = [control_output_path]

    # ✅ Remove each variable one at a time
    for removed_file in files_list_dic:
        remaining_files = [f for f in files_list_dic if f != removed_file]
        remov = removed_file

        # ✅ Reweight remaining files so that they sum to 1
        num_remaining = len(remaining_files)
        for file_info in remaining_files:
            file_info['Weight'] = 1 / num_remaining
            if num_remaining == 0:
                num_remaining = 1e-6  # Add small epsilon to avoid division by zero

        output_filename = f"{run_name}_RM_{removed_file['File'].replace('.tif', '')}.tif"
        print(f"Generating scenario without {removed_file['File']} -> {output_filename}")

        output_path = create_weighted_index(remaining_files, output_filename, sensitivity_dir)

        if output_path is None or not os.path.exists(output_path):
            print(f"❌ ERROR: Expected file {output_filename} was not created. Skipping.")
            continue

        output_files.append(output_path)

        plot_histogram(output_path,remov)  # ✅ Plot histogram for removed-variable scenario

        # ✅ Read and analyze output raster
        with rasterio.open(output_path) as src:
            data = src.read(1)
            data = data[~np.isnan(data)]
            if data.size == 0:
                print(f"❌ Warning: No valid data in {output_filename}. Skipping.")
                continue

            percent_above_5 = np.mean(data > 5) * 100
            percent_above_6 = np.mean(data > 6) * 100
            percent_above_7 = np.mean(data > 7) * 100
            d = data
            d[d == -9999] = np.nan
            nanmean_value = np.nanmean(d)

        results.append([removed_file["File"], percent_above_5, percent_above_6, percent_above_7, nanmean_value])

    # ✅ Create DataFrame and Save Results
    columns = ["Removed Variable", "% > 5", "% > 6", "% > 7", "Mean"]
    results_df = pd.DataFrame(results, columns=columns)

    print("\nSensitivity Analysis Results:")
    print(results_df)

    # Save as a single combined CSV file
    results_csv_path = os.path.join(sensitivity_dir, f"{run_name}_RM_combined_results.csv")
    results_df.to_csv(results_csv_path, index=False)

    print(f"\n✅ Sensitivity analysis completed. Results saved at: {results_csv_path}")
    return output_files, sensitivity_dir, run_name
#-------------------------------------------------------
def generate_summary_markdown(mean_raster,std_raster,range_raster,cv_raster,sensitivity_dir,run_name, mode = "w"):
    """
    Generates a markdown table summarizing the min, max, and mean values of the summary rasters.
    """
    # ✅ Mask NoData values (-9999) or other large negative placeholders
    mean_raster1 = np.where(mean_raster == -9999, np.nan, mean_raster)
    std_raster1 = np.where(std_raster == -9999, np.nan, std_raster)
    range_raster1 = np.where(range_raster == -9999, np.nan, range_raster)
    cv_raster1 = np.where(cv_raster == -9999, np.nan, cv_raster)
    stats = {
        'Raster': ['Mean', 'Standard Deviation', 'Range', 'Coefficient of Variation'],
        'Min': [
            np.nanmin(mean_raster1),
            np.nanmin(std_raster1),
            np.nanmin(range_raster1),
            np.nanmin(cv_raster1)
        ],
        'Max': [
            np.nanmax(mean_raster1),
            np.nanmax(std_raster1),
            np.nanmax(range_raster1),
            np.nanmax(cv_raster1)
        ],
        'Mean': [
            np.nanmean(mean_raster1),
            np.nanmean(std_raster1),
            np.nanmean(range_raster1),
            np.nanmean(cv_raster1)
        ]
    }

    # Convert to DataFrame
    df = pd.DataFrame(stats)

    # Format as markdown
    markdown_table = df.to_markdown(index=False, tablefmt="pipe")

    match mode:
        case 'w':
            # Define output file path
            output_path = os.path.join(sensitivity_dir, f"{run_name}_w_summary_stats.md")
        case 'r':
            # Define output file path
            output_path = os.path.join(sensitivity_dir, f"{run_name}_r_summary_stats.md")

    # Write to file
    with open(output_path, 'w') as f:
        f.write(f"## Summary Statistics for {mode}-{run_name}\n\n")
        f.write(markdown_table)
        f.write("\n\n")

    print(f"✅ Summary statistics markdown table saved at: {output_path}")
#-------------------------------------------------------
def generate_summary_rasters(output_files, sensitivity_dir, run_name, mode = 'w'):
    """
        Generates summary rasters (mean, standard deviation, range, CV) from a list of input rasters.

        Args:
            output_files (list): List of file paths to raster files.
            sensitivity_dir (str): Path to the directory where output rasters will be saved.
            run_name (str): Name for the run (used in output file names).
            mode (str, optional): sensitivity analysis mode of input ('w' for weights, 'r' for removal). Defaults to 'w'.

        Returns:
            None
        """
    # Read all rasters and stack into 3D array
    stack = []
    profile = None

    for path in output_files:
        with rasterio.open(path) as src:
            if profile is None:
                profile = src.profile.copy()
            data = src.read(1).astype(np.float32)

            # ✅ Mask NoData values (but keep it as NaN, don't zero it out yet)
            nodata_value = src.nodata
            if nodata_value is not None:
                data = np.where(data == nodata_value, np.nan, data)

            stack.append(data)

    stack = np.stack(stack)

    # ✅ Check if stack is empty before proceeding
    if stack.size == 0 or np.all(np.isnan(stack)):
        raise ValueError("Stack is empty or contains only NaN values.")

    # ✅ Check for consistent raster shapes
    for i, data in enumerate(stack):
        if data.shape != stack[0].shape:
            raise ValueError(f"Raster at index {i} has mismatched shape: {data.shape}")

    # ✅ Calculate summary statistics
    mean_raster = np.nanmean(stack, axis=0)
    std_raster = np.nanstd(stack, axis=0)
    range_raster = np.nanmax(stack, axis=0) - np.nanmin(stack, axis=0)

    # ✅ Handle division by zero in CV calculation
    with np.errstate(divide='ignore', invalid='ignore'):
        cv_raster = np.divide(std_raster, mean_raster, where=(mean_raster != 0))
        cv_raster[~np.isfinite(cv_raster)] = -9999  # Replace NaN and inf with -9999

    ## ✅ Mask the results where mean_raster is NaN (AFTER computing stats)
    #mean_raster = np.where(np.isnan(mean_raster), -9999, mean_raster)
    #std_raster = np.where(np.isnan(std_raster), -9999, std_raster)
    #range_raster = np.where(np.isnan(range_raster), -9999, range_raster)
    #cv_raster = np.where(np.isnan(cv_raster), -9999, cv_raster)

    # ✅ Function to save rasters
    def save_raster(data, name, profile):
        output_path = os.path.join(sensitivity_dir, f"{run_name}_{name}.tif")
        profile.update(dtype=rasterio.float32, nodata=-9999)
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(data, 1)
        print(f"✅ {name} raster saved at: {output_path}")

    # ✅ Set output filename prefix
    prf = "w" if mode == 'w' else "r"

    # ✅ Save rasters
    save_raster(mean_raster, f"{prf}_mean", profile)
    save_raster(std_raster, f"{prf}_std", profile)
    save_raster(range_raster, f"{prf}_range", profile)
    save_raster(cv_raster, f"{prf}_cv", profile)

    # ✅ Generate markdown summary after fixing NoData handling
    generate_summary_markdown(mean_raster, std_raster, range_raster, cv_raster, sensitivity_dir, run_name, prf)
#-----------------------------------------------------------
def generate_markdown_from_csv(csv_filepath, index_name):
    """
    Generate markdown tables from a combined results CSV file.

    Args:
        csv_filepath (str): Path to the combined results CSV file.
        index_name (str): Name of the sensitivity index (e.g., VDI, SMII, SbII, Final).

    Returns:
        str: Path to the generated markdown file.
    """
    # ✅ Read the CSV into a DataFrame
    df = pd.read_csv(csv_filepath)

    # ✅ Define output markdown path
    output_path = csv_filepath.replace('.csv', '_tables.md')

    # ✅ Create the markdown table for weight permutations
    header = f"## Sensitivity analysis on {index_name} where variable weights were swapped\n\n"
    table = df.to_markdown(index=False, tablefmt="pipe")

    # ✅ Bold the first row (assuming it's the primary analysis)
    table = table.replace(df.iloc[0, 0], f"**{df.iloc[0, 0]}**")
    for col in df.columns[1:]:
        table = table.replace(str(df.iloc[0][col]), f"**{df.iloc[0][col]}**")

    # ✅ Add a descriptive caption
    caption = f"\n: Sensitivity analysis on {index_name} with swapped variable weights. " \
              f"The values from the row in bold were used in the primary analysis.\n"

    # ✅ Create the markdown table for removal scenarios if applicable
    if 'Removed Variable' in df.columns:
        removal_header = f"\n## Sensitivity analysis on {index_name} with removal of variables (OAT)\n\n"
        removal_table = df.to_markdown(index=False, tablefmt="pipe")

        # Bold the "None" row (assuming it's the baseline)
        removal_table = removal_table.replace('None', '**None**')
        for col in df.columns[1:]:
            baseline_row = df.loc[df['Removed Variable'].str.lower() == 'none']
            if not baseline_row.empty:
                for c in df.columns[1:]:
                    removal_table = removal_table.replace(
                        str(baseline_row[c].values[0]),
                        f"**{baseline_row[c].values[0]}**"
                    )

        removal_caption = f"\n: Sensitivity analysis on {index_name} with removal of variables One-at-a-Time (OAT). " \
                          f"The values from the row in bold were used in the primary analysis.\n"
    else:
        removal_header = ""
        removal_table = ""
        removal_caption = ""

    # ✅ Write to markdown file
    with open(output_path, 'w') as f:
        f.write(header)
        f.write(table)
        f.write(caption)
        if removal_header:
            f.write(removal_header)
            f.write(removal_table)
            f.write(removal_caption)

    print(f"✅ Markdown table saved at: {output_path}")
    return output_path
####################### testing --------------------------------------------------

print("SuitabilityMappingTools has been loaded")


# test generate_markdown_from_csv()
#path = r'Data/Sensitivity_Analysis/SbII/SbII_RM_combined_results.csv'
#generate_markdown_from_csv(path,'SbII')