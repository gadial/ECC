/* 
 * File:   hcp.h
 * Author: gadial
 *
 * Created on December 3, 2009, 11:50 AM
 */

#ifndef _HCP_H
#define	_HCP_H

#include <map>
#include <string>
using std::map;
using std::string;
        
class HCP{ //Hilbert class polynomials
public:
    HCP();
private:
    map<int, string> H;
};


#endif	/* _HCP_H */

