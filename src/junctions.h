/*
  Anything related to junctions are goes in here
*/

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <any>
#include <set>

struct Junctions;
struct SingleJunction;

int 
take_closest (std::vector<int> list, int num);

Junctions 
blocks_to_junctions (std::vector<std::pair<int, int>> blocks);

Junctions 
get_TSS_TES_site (std::map<std::string, Junctions> transcript_to_junctions, std::vector<std::string> tr_list);

std::set<int> 
get_splice_site (std::map<std::string, Junctions> transcript_to_junctions, std::vector<std::string> tr_list);