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
from rasterio.plot import plotting_extent
import pandas as pd
from rasterio.warp import reproject, Resampling
import rasterio
from rasterio.plot import show, plotting_extent
from rasterio.enums import Resampling
import os
import glob
from rasterio.plot import show
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import seaborn as sns

# Functions---------------------------------
def plot_layers(data_dir):
    # Get list of all raster files in the data directory
    raster_files = [f for f in os.listdir(data_dir) if f.endswith('.tif')]

    # Set up the plot
    num_files = len(raster_files)
    fig, axes = plt.subplots(1, num_files, figsize=(15, 5))

    if num_files == 1:
        axes = [axes]

    for ax, raster_file in zip(axes, raster_files):
        with rasterio.open(os.path.join(data_dir, raster_file)) as src:
            data = src.read(1)
            # Mask no-data values
            data = data.astype(float)
            data[data == src.nodata] = float('nan')

            # Plot the data with scaling based on min and max values
            img = ax.imshow(data, cmap='BrBG', vmin=np.nanmin(data), vmax=np.nanmax(data))
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
def plot_histogram(filepath, weights):
    """
    This function loads a raster file, plots a histogram of values with 10 bins between 1 and 10,
    provides a count of pixels excluding no-data pixels, and calculates the percentage of pixels
    in each bin. The figure is saved as a .png file with the same basename as the filepath provided.
    Additionally, it saves a markdown table as a .txt file showing the percentage of total pixels
    (excluding no-data) in each bin, the count of pixels in each bin, and the min and max values of each bin.
    """
    with rasterio.open(filepath) as src:
        data = src.read(1)
        nodata = src.nodata

        # Exclude no-data pixels
        if nodata is not None:
            data = data[data != nodata]

        # Define bins
        bins = np.arange(0, 11, 1)

        # Plot histogram
        fig, ax = plt.subplots()
        counts, bins, patches = ax.hist(data, bins=bins, edgecolor='black')


        ax.set_title(f"Histogram of {os.path.basename(filepath)}: {weights}")
        ax.set_xlabel("Value")
        ax.set_ylabel("Frequency (thousands of pixels)")
        ax.set_xticks(bins)
        ax.set_yticklabels([f'{int(label/1000)}k' for label in ax.get_yticks()])

        # Count of pixels excluding no-data pixels
        pixel_count = np.sum(data > 0.5)

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
            f.write("| Bin | Count | Percentage of Total Pixels | Min Value | Max Value | Area ha |\n")
            f.write("|-----|-------|----------------------------|-----------|-----------|---------|\n")
            for i in range(len(counts)):
                bin_min = bins[i]
                bin_max = bins[i + 1]
                pixels = int(counts[i])
                sqmet = pixels * 900
                hectares = sqmet/10000
                f.write(f"| {int(bins[i])} | {int(counts[i])} | {percentages[int(bins[i])]:.2f}% | {bin_min:.2f} | {bin_max:.2f} |{hectares}|\n")

        print(f"Count of pixels (excluding no-data pixels): {pixel_count}")
        for i in range(len(counts)):
            print(f"Percentage of pixels in bin {int(bins[i])}: count = {int(counts[i])} :{percentages[int(bins[i])]:.2f}%")
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

    weighted_sum[weighted_sum < 0.001] = np.nan

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

    sensitivity_dir = os.path.join(os.path.dirname(files_list_dic[0]['Filepath']), "Sensitivity_Analysis")
    os.makedirs(sensitivity_dir, exist_ok=True)

    file_names = [file_info['File'] for file_info in files_list_dic]
    base_weights = [file_info['Weight'] for file_info in files_list_dic]
    weight_permutations = list(itertools.permutations(base_weights))

    results = []

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

        # Analyze the newly created raster
        with rasterio.open(output_path) as src:
            data = src.read(1)
            nodata = src.nodata

        if nodata is not None:
            data = data[data != nodata]  # Remove NoData values

        percent_above_5 = np.mean(data > 5) * 100
        percent_above_6 = np.mean(data > 6) * 100
        percent_above_7 = np.mean(data > 7) * 100
        d = data
        d[d == -9999] = np.nan
        nanmean_value = np.nanmean(d)

        results.append([output_filename] + list(weight_set) + [percent_above_5, percent_above_6, percent_above_7, nanmean_value])

    # Save results
    results_df = pd.DataFrame(results, columns=["File", *file_names, "% > 5", "% > 6", "% > 7", "Mean"])
    results_txt_path = os.path.join(sensitivity_dir, f"{run_name}_WSA_results.txt")
    results_df.to_csv(results_txt_path, sep="\t", index=False)

    print(f"✅ Sensitivity analysis completed. Results saved at: {results_txt_path}")
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
    sensitivity_dir = os.path.join(os.path.dirname(files_list_dic[0]['Filepath']), "Sensitivity_Analysis")
    os.makedirs(sensitivity_dir, exist_ok=True)

    # ✅ Ensure define_weights() returns a valid list
    control_weights = define_weights(files_list_dic, False)

    if control_weights is None or not isinstance(control_weights, list) or len(control_weights) == 0:
        print(f"❌ ERROR: `define_weights()` returned an invalid list. Check input files.")
        return

    # ✅ Control scenario with all variables equally weighted
    control_output_filename = f"{run_name}_RM_control.tif"
    print(f"Generating control scenario: {control_output_filename}")

    control_output_path = create_weighted_index(control_weights, control_output_filename, sensitivity_dir)

    if control_output_path is None or not os.path.exists(control_output_path):
        print(f"❌ ERROR: Expected control file {control_output_filename} was not created. Skipping analysis.")
        return

    plot_histogram(filepath=control_output_path,weights= remov)  # ✅ Plot histogram for control scenario

    # ✅ Read and analyze control raster
    with rasterio.open(control_output_path) as src:
        data = src.read(1)
        data = data[~np.isnan(data)]

        control_percent_above_5 = np.mean(data > 5) * 100
        control_percent_above_6 = np.mean(data > 6) * 100
        control_percent_above_7 = np.mean(data > 7) * 100
        d = data
        d[d == -9999] = np.nan
        control_nanmean_value = np.nanmean(d)

    results = [
        ["Control", control_percent_above_5, control_percent_above_6, control_percent_above_7, control_nanmean_value]
    ]

    # ✅ Remove each variable one at a time
    for removed_file in files_list_dic:
        remaining_files = [f for f in files_list_dic if f != removed_file]
        remov = removed_file

        # ✅ Ensure valid weights after removal
        adjusted_weights = define_weights(remaining_files)

        if adjusted_weights is None or len(adjusted_weights) == 0:
            print(f"❌ ERROR: No valid weights after removing {removed_file['File']}. Skipping.")
            continue

        output_filename = f"{run_name}_RM_{removed_file['File'].replace('.tif', '')}.tif"
        print(f"Generating scenario without {removed_file['File']} -> {output_filename}")

        output_path = create_weighted_index(adjusted_weights, output_filename, sensitivity_dir)

        if output_path is None or not os.path.exists(output_path):
            print(f"❌ ERROR: Expected file {output_filename} was not created. Skipping.")
            continue

        plot_histogram(output_path,remov)  # ✅ Plot histogram for removed-variable scenario

        # ✅ Read and analyze output raster
        with rasterio.open(output_path) as src:
            data = src.read(1)
            data = data[~np.isnan(data)]

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

    results_txt_path = os.path.join(sensitivity_dir, f"{run_name}_RM_results.txt")
    results_df.to_csv(results_txt_path, sep="\t", index=False)

    print(f"\n✅ Sensitivity analysis completed. Results saved at: {results_txt_path}")
####################### testing --------------------------------------------------

print("SuitabilityMappingTools has been loaded")
