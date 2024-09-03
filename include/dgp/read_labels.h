#ifndef LIBDGP_READ_LABELS_H
#define LIBDGP_READ_LABELS_H

#include "dgp_inline.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

namespace dgp
{
    template<typename T> 
    DGP_INLINE bool read_labels(std::string filename, std::vector<T> &labels, int skip_lines = 1)
    {
        std::ifstream ifs(filename);
        if (!ifs.is_open())
        {
            return false;
        }

        labels.clear();
        std::string line;
        for(;std::getline(ifs, line);)
        {
            if(skip_lines-- > 0)
            {
                continue;
            }

            std::istringstream iss(line);
            T label;
            iss >> label;
            labels.push_back(label);
        }

        ifs.close();
        return true;
    }
}

#endif //LIBDGP_READ_LABELS_H