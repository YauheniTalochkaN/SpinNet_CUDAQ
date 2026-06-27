#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

inline void export_csv_helper(std::vector<std::string>& headers,
                              std::vector<const std::vector<double>*>& columns) 
{

}

template <typename... Rest>
inline void export_csv_helper(std::vector<std::string>& headers,
                              std::vector<const std::vector<double>*>& columns,
                              const std::pair<std::string, const std::vector<double>&>& data,
                              const Rest&... rest) 
{
    headers.push_back(data.first);
    columns.push_back(&data.second);

    export_csv_helper(headers, columns, rest...);
}

template <typename... Args>
void export_csv(const std::string& filename, 
                const std::pair<std::string, const std::vector<double>&>& fdata,
                const Args&... args) 
{
    std::vector<std::string> headers;
    std::vector<const std::vector<double>*> columns;

    headers.push_back(fdata.first);
    columns.push_back(&fdata.second);

    export_csv_helper(headers, columns, args...);

    size_t n = columns.front()->size();

    for(const auto* column : columns) 
    {
        if(column->size() != n) 
        {
            std::cerr << "Error: all columns must have the same length." << std::endl;
            return;
        }
    }

    std::ofstream file(filename);

    if(!file.is_open()) 
    {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    for(size_t i = 0; i < headers.size(); ++i) 
    {
        file << headers[i];

        if(i < headers.size() - 1) 
        {
            file << ",";
        }
    }

    file << std::endl;

    for(size_t i = 0; i < n; ++i) 
    {
        for(size_t j = 0; j < columns.size(); ++j) 
        {
            file << std::fixed << std::setprecision(8) << (*columns[j])[i];
            
            if(j < columns.size() - 1) 
            {
                file << ",";
            }
        }

        file << std::endl;
    }
    
    file.close();
}