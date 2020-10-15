//
// Created by Tobi Popoola on 10/6/20.
//

#include "Utils.h"

void Utils::replaceAllString(std::string &Original,
                             const std::string &ToReplace,
                             const std::string &Replacement) {
    size_t pos = Original.find(ToReplace);
    while (pos != std::string::npos) {
        Original.replace(pos, ToReplace.length(), Replacement);
        pos = Original.find(ToReplace);
    }
}
