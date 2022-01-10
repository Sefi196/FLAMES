/*
  miscellaneous helper functions and things
  they will not live here forever, 
  once I figure out how to categorize them better then I will move them
*/
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "../classes/Pos.h"
#include "../classes/StartEndPair.h"

#ifndef ISO
#define ISO

struct Iso
{
  /*
    a data container used in Isoforms,
    specifically for known_isoforms and match_known_annotation
  */

  long support_count;
  std::string transcript_id; 
  std::string gene_id;
};

#endif

std::vector<StartEndPair>
pairwise (std::vector<int> input);

std::vector<int>
find_best_splice_chain(std::vector<int> raw_iso, std::vector<std::vector<int>> junction_list, int MAX_DIST);

int
if_exon_contains(std::vector<int> s1, std::vector<int> s2, int MAX_TOLERANCE);

float
get_exon_sim_pct(std::vector<int> exons1, std::vector<int> exons2);

int
iv_overlap (StartEndPair iv1, StartEndPair iv2);

int
exon_overlap (std::vector<int> exons1, std::vector<StartEndPair> exons2);

std::vector<std::pair<std::string, std::string>>
get_fa(std::string filename);