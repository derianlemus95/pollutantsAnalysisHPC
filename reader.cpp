#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <algorithm>

namespace fs = std::filesystem;
int verbose = 0;

std::vector<std::vector<std::string>> readCSV(const fs::path& filepath) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filepath);

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            row.push_back(value);
        }
        if(row[4].find("-999.0") == std::string::npos){
            data.push_back(row);
        } 
    }
    return data;
}

float convertToUGM3(float value, const std::string& unit, const std::string& substance) {
    float molecularWeight = 0.0;
    if (substance == "NO2") molecularWeight = 46.0055;
    else if (substance == "SO2") molecularWeight = 64.066;
    else if (substance == "CO") molecularWeight = 28.01;
    else if (substance == "OZONE") molecularWeight = 48.00;
    // PM2.5 and PM10 are typically measured directly in µg/m³, so no conversion is needed.

    if (unit == "ppm") {
        return value * molecularWeight * 1000000 / 24.45;
        //return value * molecularWeight * 1000 / 24.45;
    } else if (unit == "ppb") {
        return value * molecularWeight * 1000 / 24.45;
        //return value * molecularWeight  / 24.45;
    }
    return value; // Assume the value is already in µg/m³ or no conversion needed
}

void writeToCSV(const fs::path& filepath, std::vector<std::vector<std::string>>& data, const std::vector<float>& zScores) {
    //Create Zscores column and add at the back of the columns
    //for(size_t i =0; i < data.size(); ++i){
      //  data[i].push_back(zScores[i]);
   // }
    
    //for(auto& row: data){
      //  std::cout << row[13] << " " ;
    //}   
             
    //}
    //modified wrting to file since it was causing errors
    //now only writes the substance and concentration to the ouput file 
    std::ofstream file(filepath, std::ios::trunc);
    for (size_t i = 0; i < data.size(); ++i) {
            if(data[i].size() > 3) {
                file << data[i][3] << "," << std::to_string(zScores[i]);
                file << "\n";
            }	
    }
        
    
    file.close();
}

void normalizeAndConvertFile(const fs::path& inputPath, const fs::path& outputPath) {
    auto data = readCSV(inputPath);
    std::vector<float> concentrations;
    for (size_t i = 0; i < data.size(); ++i) { // Skip header row
        try {
            std::string x = data[i][4];
            x.erase(std::remove(x.begin(), x.end(), '"'), x.end());
            //std::cout << x;
            float value = std::stof(x); // store concentration as a float
            //if value is invalid store -999 for concentration, else calculate concentration
            if(value != -999){
                //std::cout << value << " " ;
                const std::string& unit = data[i][5]; 
                const std::string& substance = data[i][3]; 
                concentrations.push_back(convertToUGM3(value, unit, substance));
            } 
            // else {
            //     const std::string& unit = data[i][5]; 
            //     const std::string& substance = data[i][3]; 
            //     concentrations.push_back(-999.0);
            // }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Conversion issue on line " << i << " (" << e.what() << ")" << std::endl;
            // Insert a default or error value if conversion fails
            concentrations.push_back(0.0);
        }
    }

    //filter out invalid concentrations
    auto tempVector = concentrations;
    
    tempVector.erase(std::remove(tempVector.begin(), tempVector.end(), -999.0), tempVector.end());
    
    // Calculate the mean and standard deviation for Z-Score normalization
    float mean = std::accumulate(tempVector.begin(), tempVector.end(), 0.0) / tempVector.size();
    std::vector<float> diff(tempVector.size());
    std::transform(tempVector.begin(), tempVector.end(), diff.begin(), [mean](float x) { return x - mean; });
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / tempVector.size());
 
    //store z-score normalizatoion in a vector
    std::vector<float> zScores;
    // Apply Z-Score normalization to the concentration values
    for (size_t i = 0; i < concentrations.size(); ++i) {
        if(concentrations[i] != -999.0){
            float zScore = (concentrations[i] - mean) / stdev;
            //zScores.push_back(std::to_string(zScore));
            zScores.push_back(zScore);
        } 
        // else {
        //     zScores.push_back(std::to_string(-999.0));
        // }
    }
    //for(auto& row: zScores){
   //     if(row > 10){
   //     std::cout << row << " " ;
   // }
    //}
    writeToCSV(outputPath, data, zScores);
}

void findCSVFiles(const fs::path& rootPath, std::vector<fs::path>& outFilePaths) {
    for (const auto& entry : fs::recursive_directory_iterator(rootPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv" && entry.path().string().find("norm") == std::string::npos) {
            outFilePaths.push_back(entry.path());
        }
    }
}

int main(int argc, char** argv) {
    //CLANG: std::__1::chrono::steady_clock::time_point start;
    //GCC
    std::chrono::steady_clock::time_point start;
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::vector<fs::path> allCSVFiles;
    if (world_rank == 0) {
        std::cout << "Starting Parallel File processing...." << std::endl;
        start = std::chrono::steady_clock::now();
        // Master process: Find CSV files and distribute file paths to workers
        findCSVFiles("./data", allCSVFiles);
        
        int num_files = allCSVFiles.size();
        for (int i = 0; i < num_files; ++i) {
            int target_process = i % (world_size - 1) + 1; // Distribute files to workers, starting from process 1
            std::string filepath = allCSVFiles[i].string();
            MPI_Send(filepath.c_str(), filepath.size() + 1, MPI_CHAR, target_process, 0, MPI_COMM_WORLD);
        }

        // Sending an empty string as a signal to stop
        for (int i = 1; i < world_size; ++i) {
            char end[3] = {'E','N','D'};
            if(verbose > 0)
                std::cout << "Sending breaker to: " << i << std::endl;
            MPI_Send(end, 3, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // Worker processes: Receive file paths and process them
        try {
            auto count = 0;
            char filepath[256];
            MPI_Status status;
            while (true) {
                // Receive the actual message
                MPI_Recv(filepath, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (std::string(filepath) == "END") {
                    if(verbose > 0)
                        std::cout << "Rank: " << world_rank << ", Received breaker" << ", Stopping..." << std::endl;
                    break;
                } // Termination signal received
                count++;
                if(verbose > 0)
                    std::cout << "Rank: " << world_rank << ", received file for processing: " << filepath << ", total count: " << count << std::endl;
                fs::path inputPath(filepath);
                fs::path outputPath = inputPath.parent_path() / ("normalized_" + inputPath.filename().string());
                normalizeAndConvertFile(inputPath, outputPath);
                memset(filepath, 0, 256);
            }   
        } catch(const std::exception& e) {
            std::cerr << "Exception" << e.what() << std::endl;
        }
        
    }
    MPI_Finalize();
    auto finish = std::chrono::steady_clock::now();
    if(world_rank == 0)
        std::cout << "Processed: " << allCSVFiles.size() << " files in: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() / 1000000000.0 << "s" << std::endl;
    return 0;
}
