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

/*Structure:
	1. Master rank distributes the work to other processes
	2. Each worker process reads a normalized.csv file
	1. As its reading it, it filters by substance and stores its corresponding norm value to its local substance vector
	3. Return these local vectors for each substance and sends them back to master rank
	4. Master rank recieves each vector from each process and appends it to the global vector for each substance
	5. Calculate the min size of the global vectors
	6. Use min size for computation of correlation for each substance
	7. Output the results to a csv file
	   a. Contents are in the form 
	          PM2.5  PM10  OZONE  NO2  CO SO2
	   PM2.5
	   PM10
	   OZONE
	   NO2
	   CO
	   SO2
*/

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
        if(row[1].find("-999.0") == std::string::npos){
            data.push_back(row);
        } 
    }
  
  //for(auto& row: data){
    // for(auto& x: row){
      //       std::cout << x << " " ;
     //}    
             
    //}
    //return only substance and normalized columns
    return data;
}

//function to write correlation matrix onto csv file
void writeToCSV(const std::vector<float>& v1, const std::vector<float>& v2, const std::vector<float>& v3,    
const std::vector<float>& v4, const std::vector<float>& v5, const std::vector<float>& v6){

const std::string& fileName = "correlationMatrix.csv";
std::ofstream file(fileName, std::ios::out | std::ios::trunc);

for (size_t i = 0; i < v1.size(); i++){
    file << std::to_string(v1[i])<< "," << std::to_string(v2[i]) << "," << std::to_string(v3[i]) << "," << std::to_string(v4[i]) << "," << std::to_string(v5[i]) << "," << std::to_string(v6[i]);
    file << "\n";
}
file.close();
}


float computeCorrelation(std::vector<float> s1, std::vector<float> s2, int n){
float sum1 = 0.0, sum2 = 0.0, sum1Times2 = 0.0, sum1Squared=0.0, sum2Squared=0.0;
for (int i =0; i < n; i++){
    sum1 += s1[i];
    sum2 += s2[i];
    sum1Times2 +=  s1[i] * s2[i];
    sum1Squared += s1[i] * s1[i];
    sum2Squared += s2[i] * s2[i];    
}

float numerator = n * sum1Times2 - sum1 * sum2;
float denominator = std::sqrt((n* sum1Squared - sum1 * sum1) * (n * sum2Squared - sum2 * sum2));

return numerator/ denominator;
}

std::vector<float> filterBySubstance(const fs::path& inputPath, std::string substance) {
    auto data = readCSV(inputPath);
    //for(auto& row: data){
       // std::cout << row[1] << " " ;
   //}
    std::vector<float> tempData;
    for (size_t i = 0; i < data.size(); ++i) { 
        try {               
               //data[i][3] is different type than "PM2.5"
                std::string& s = data[i][0];
                s.erase(std::remove(s.begin(), s.end(), '"'), s.end());
                //std::cout << s << " (" << substance << " )";
                if(s == substance){
                
                std::string& z = data[i][1];
                z.erase(std::remove(z.begin(), z.end(), '"'), z.end());
                
                float normalized = std::stof(z);
                //round norm value to 3 spaces after decimal
                //float roundedNumber = std::round(normalized * 1000)/ 1000;
                //std::cout << roundedNumber << " " ;
                tempData.push_back(normalized);
  
            }      
        } catch (const std::exception& e) {
            std::cerr << "Warning: Filtering data error on line " << i << " (" << e.what() << ")" << std::endl;
        }      
    }
    return tempData;
}



void findCSVFiles(const fs::path& rootPath, std::vector<fs::path>& outFilePaths) {
    for (const auto& entry : fs::recursive_directory_iterator(rootPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv" && entry.path().filename().string().rfind("normalized_", 0) == 0) {
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
    
    std::vector<float> globalPm2Data;
    std::vector<float> globalPm10Data;
    std::vector<float> globalOzoneData;
    std::vector<float> globalNo2Data;
    std::vector<float> globalCoData;
    std::vector<float> globalSo2Data;
    
    if (world_rank == 0) {
        std::cout << "Starting Parallel File processing...." << std::endl;
        start = std::chrono::steady_clock::now();
        // Master process: Find CSV files and distribute file paths to workers
        findCSVFiles("./data", allCSVFiles);
        //findCSVFiles("./data", allCSVFiles);
        
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
        //recieve localPm2Vectors and append them to the global Vector
        for(int i = 1; i< world_size; i++){
            int dataSize;
            MPI_Recv(&dataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<float> tempData(dataSize);
            MPI_Recv(tempData.data(), dataSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            globalPm2Data.insert(globalPm2Data.end(), tempData.begin(), tempData.end());
        }
        //recieve localPm10Vectors and append them to the global Vector
        for(int i = 1; i< world_size; i++){
            int dataSize;
            MPI_Recv(&dataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<float> tempData(dataSize);
            MPI_Recv(tempData.data(), dataSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            globalPm10Data.insert(globalPm10Data.end(), tempData.begin(), tempData.end());
        }
        //recieve localOzoneVectors and append them to the global Vector
        for(int i = 1; i< world_size; i++){
            int dataSize;
            MPI_Recv(&dataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<float> tempData(dataSize);
            MPI_Recv(tempData.data(), dataSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            globalOzoneData.insert(globalOzoneData.end(), tempData.begin(), tempData.end());
        }
        //recieve localNo2Vectors and append them to the global Vector
        for(int i = 1; i< world_size; i++){
            int dataSize;
            MPI_Recv(&dataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<float> tempData(dataSize);
            MPI_Recv(tempData.data(), dataSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            globalNo2Data.insert(globalNo2Data.end(), tempData.begin(), tempData.end());
        }
        //recieve localCoVectors and append them to the global Vector
        for(int i = 1; i< world_size; i++){
            int dataSize;
            MPI_Recv(&dataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<float> tempData(dataSize);
            MPI_Recv(tempData.data(), dataSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            globalCoData.insert(globalCoData.end(), tempData.begin(), tempData.end());
        }
        //recieve localSo2Vectors and append them to the global Vector
        for(int i = 1; i< world_size; i++){
            int dataSize;
            MPI_Recv(&dataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<float> tempData(dataSize);
            MPI_Recv(tempData.data(), dataSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            globalSo2Data.insert(globalSo2Data.end(), tempData.begin(), tempData.end());
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
               	memset(filepath, 0, 256); 
               	
               	//std::cout << localPm2Data.size() << std::endl;
                          	
               	//send back localPm2Data 
               	std::vector<float> localPm2Data = filterBySubstance(inputPath, "PM2.5");
               	int localPm2DataSize = localPm2Data.size();//chat gtp
               	MPI_Send(&localPm2DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);//chatgpt
               	MPI_Send(localPm2Data.data(), localPm2DataSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);//chatgpt
                //send back localPm10Data 
               	std::vector<float> localPm10Data = filterBySubstance(inputPath, "PM10");
               	int localPm10DataSize = localPm10Data.size();
               	MPI_Send(&localPm10DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
               	MPI_Send(localPm10Data.data(), localPm10DataSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);	
                
                //send back localOzoneData 
                std::vector<float> localOzoneData = filterBySubstance(inputPath, "OZONE");
                int localOzoneDataSize = localOzoneData.size();
               	MPI_Send(&localOzoneDataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
               	MPI_Send(localOzoneData.data(), localOzoneDataSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);	
                
                //send back localNo2Data
                std::vector<float> localNo2Data = filterBySubstance(inputPath, "NO2");
                int localNo2DataSize = localNo2Data.size();
               	MPI_Send(&localNo2DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
               	MPI_Send(localNo2Data.data(), localNo2DataSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);	
                              
                //send back localCoData 
                std::vector<float> localCoData = filterBySubstance(inputPath, "CO");
                int localCoDataSize = localCoData.size();
               	MPI_Send(&localCoDataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
               	MPI_Send(localCoData.data(), localCoDataSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);	
                              
                //send back localSo2Data 
                std::vector<float> localSo2Data = filterBySubstance(inputPath, "SO2");
                int localSo2DataSize = localSo2Data.size();
               	MPI_Send(&localSo2DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
               	MPI_Send(localSo2Data.data(), localSo2DataSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); 	
     
            }   
        } catch(const std::exception& e) {
            std::cerr << "Exception" << e.what() << std::endl;
        }
        
    }
    
    MPI_Finalize();
    auto finish = std::chrono::steady_clock::now();
    if(world_rank == 0){
        std::cout << "Processed: " << allCSVFiles.size() << " files in: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() / 1000000000.0 << "s" << std::endl;
        std::cout << "ALl PM2Data: " << globalPm2Data.size() << std::endl;
        std::cout << "ALl PM10Data: " << globalPm10Data.size() << std::endl;
        std::cout << "ALl OzoneData: " << globalOzoneData.size() << std::endl;
        std::cout << "ALl No2Data: " << globalNo2Data.size() << std::endl;
        std::cout << "ALl CoData: " << globalCoData.size() << std::endl;
        std::cout << "ALl So2Data: " << globalSo2Data.size() << std::endl;
         
         //get lowest size value
         size_t n = std::min({globalPm2Data.size(), globalPm10Data.size(), globalOzoneData.size(), globalNo2Data.size(), globalCoData.size(), globalSo2Data.size()});
         std::cout << n << std::endl;
         
         std::vector<float> pm2Row;
         pm2Row.push_back(1.00); 
         pm2Row.push_back(computeCorrelation( globalPm2Data, globalPm10Data, n));
         pm2Row.push_back(computeCorrelation( globalPm2Data, globalOzoneData, n));
         pm2Row.push_back(computeCorrelation( globalPm2Data, globalNo2Data, n));
         pm2Row.push_back(computeCorrelation( globalPm2Data, globalCoData, n));
         pm2Row.push_back(computeCorrelation( globalPm2Data, globalSo2Data, n));
         
         std::vector<float> pm10Row;
         pm10Row.push_back(computeCorrelation( globalPm10Data, globalPm2Data, n));
         pm10Row.push_back(1.00); 
         pm10Row.push_back(computeCorrelation( globalPm10Data, globalOzoneData, n));
         pm10Row.push_back(computeCorrelation( globalPm10Data, globalNo2Data, n));
         pm10Row.push_back(computeCorrelation( globalPm10Data, globalCoData, n));
         pm10Row.push_back(computeCorrelation( globalPm10Data, globalSo2Data, n));
         
         std::vector<float> ozoneRow;
         ozoneRow.push_back(computeCorrelation( globalOzoneData, globalPm2Data, n));
         ozoneRow.push_back(computeCorrelation( globalOzoneData, globalPm10Data, n));
         ozoneRow.push_back(1.00); 
         ozoneRow.push_back(computeCorrelation( globalOzoneData, globalNo2Data, n));
         ozoneRow.push_back(computeCorrelation( globalOzoneData, globalCoData, n));
         ozoneRow.push_back(computeCorrelation( globalOzoneData, globalSo2Data, n));
         
         std::vector<float> no2Row;
         no2Row.push_back(computeCorrelation( globalNo2Data, globalPm2Data, n));
         no2Row.push_back(computeCorrelation( globalNo2Data, globalPm10Data, n));
         no2Row.push_back(computeCorrelation( globalNo2Data, globalOzoneData, n));
         no2Row.push_back(1.00); 
         no2Row.push_back(computeCorrelation( globalNo2Data, globalCoData, n));
         no2Row.push_back(computeCorrelation( globalNo2Data, globalSo2Data, n));
         
         std::vector<float> coRow;
         coRow.push_back(computeCorrelation( globalCoData, globalPm2Data, n));
         coRow.push_back(computeCorrelation( globalCoData, globalPm10Data, n));
         coRow.push_back(computeCorrelation( globalCoData, globalOzoneData, n));
         coRow.push_back(computeCorrelation( globalCoData, globalNo2Data, n));
         coRow.push_back(1.00); 
         coRow.push_back(computeCorrelation( globalCoData, globalSo2Data, n));
         
         std::vector<float> so2Row;
         so2Row.push_back(computeCorrelation( globalSo2Data, globalPm2Data, n));
         so2Row.push_back(computeCorrelation( globalSo2Data, globalPm10Data, n));
         so2Row.push_back(computeCorrelation( globalSo2Data, globalOzoneData, n));
         so2Row.push_back(computeCorrelation( globalSo2Data, globalNo2Data, n));
         so2Row.push_back(computeCorrelation( globalSo2Data, globalCoData, n));
         so2Row.push_back(1.00);
         
         for(auto& x: pm2Row){
             std::cout << x << " " ;
         }
         std::cout << std::endl; 
         for(auto& x: pm10Row){
             std::cout << x << " " ;
         }
         std::cout << std::endl; 
          for(auto& x: ozoneRow){
             std::cout << x << " " ;
         }
         std::cout << std::endl; 
         for(auto& x: no2Row){
             std::cout << x << " " ;
         }
         std::cout << std::endl; 
         for(auto& x: coRow){
             std::cout << x << " " ;
         }
         std::cout << std::endl; 
         for(auto& x: so2Row){
             std::cout << x << " " ;
         }
         std::cout << std::endl; 
         
         writeToCSV(pm2Row, pm10Row, ozoneRow, no2Row, coRow, so2Row);
        }
        

    return 0;
}
