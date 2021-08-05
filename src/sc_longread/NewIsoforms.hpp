/*
  the intention of NewIsoforms is to overcome problems with the original.
  Namely:
  - the single-exon/multi-exon stuff, which is wasting both space and processing
  - whatever is going on with lr_pair - what a strange choice of data structure
*/ 

#include <map>
#include <vector>
#include "junctions.hpp"

class Isoforms
{
  private:
    // these values will all be extracted from config
    const int MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
              MAX_SITE_PER_SLICE, MIN_SUP_CNT, MIN_FL_EXON_LEN,
              MIN_SUP_PCT, STRAND_SPECIFIC, REMOVE_INCOMP_READS;
  public:
    std::string ch;
    
    std::map<std::vector<int>, int> 
    junction_dict;
    
    std::vector<std::vector<std::vector<int>>> 
    junction_list;

    std::map<std::vector<int>, std::vector<int>> 
    lr_pair;
    std::vector<int> left;
    std::vector<int> right;

    std::map<std::pair<int, int>, int> single_block_dict;
    std::vector<std::vector<std::pair<int, int>>> single_blocks;
    std::map<std::pair<int, int>, std::vector<int>> strand_cnt;
    std::map<std::string, int> new_isoforms;
    std::map<std::string, int> known_isoforms;
    std::map<std::string, int> raw_isoforms;
    std::map<std::string, int> ge_dict;

    Isoforms(std::string ch, std::map<std::string, int> config);
    void add_isoform(Junctions junctions, bool is_reversed);
    void add_one(Junctions junctions, bool strand);
    void update_one(Junctions junctions, std::vector<int> key, bool strand);
};

Isoforms::Isoforms(std::string ch, std::map<std::string, int> config)
: MAX_DIST              (config["MAX_DIST"]),
  MAX_TS_DIST           (config["MAX_TS_DIST"]),
  MAX_SPLICE_MATCH_DIST (config["MAX_SPLICE_MATCH_DIST"]),
  MAX_SITE_PER_SLICE    (config["Max_site_per_splice"]),
  MIN_SUP_CNT           (config["Min_sup_cnt"]),
  MIN_FL_EXON_LEN       (config["min_fl_exon_len"]),
  MIN_SUP_PCT           (config["Min_sup_pct"]),
  STRAND_SPECIFIC       (config["strand_specific"]),
  REMOVE_INCOMP_READS   (config["remove_incomp_reads"])
{
  /*
    initialises the object,
    extracting all the const values from config
    
    (why are config keys a mix of CAPS, lower, and Sentence Case?)
    (I have no idea, that's how they were when I found them)
    (change it later)
  */
  
  this->ch = ch;
}

void Isoforms::add_isoform(Junctions junctions, bool is_reversed, bool is_single_exon=false)
{
  if (is_single_exon)
  {
    if (this->single_block_dict.size() == 0)
    {
      // this will be the first thing in single_block_dict
      // so, nothing to check
      this->add_one(junctions, is_reversed, is_single_exon);
    }
    else
    {
      // there is already stuff in single_block_dict
      // first, check if this key is already in there
      std::pair<int, int> target_key =
      {junctions.left, junctions.right};

      if (this->single_block_dict.count(target_key) > 0)
      {
        // the key is already in there

      }
    }
  }
}

void Isoforms::add_isoform_old(Junctions junctions, bool is_reversed)
{
  if (junctions.junctions.size() == 0) // single exon reads
  {
    if (this->single_block_dict.size() == 0)
    {
      // this will be the first thing in single_block_dict
      // so, nothing to check
      this->add_one(junctions, is_reversed);
    }
    else
    {
      // there is already stuff in single_block_dict
      // first, check if this key is already in there
      std::pair<int, int> target_key = 
      {junctions.left, junctions.right};

      if (this->single_block_dict.count(target_key) > 0)
      {
        // this time the key must be converted to vector
        // we really need to fix this
        std::vector<int> target_key_vector =
        {target_key.first, target_key.second};

        // the key is already present - so just update it
        this->update_one(junctions, target_key_vector, is_reversed);
      }
      else
      {
        // the key is not present in the dict yet
        // check if there's anything very close to the key
        bool found = false;

        for (auto const& [current_key, current_value] : this->single_block_dict)
        {
          if ((abs(current_key.first - target_key.first) < this->MAX_TS_DIST) and
              (abs(current_key.second - target_key.second) < this->MAX_TS_DIST))
          {
            // we've found one that's very similar
            // just update it
            this->update_one(junctions, {current_key.first, current_key.second}, is_reversed);
            found = true;
            break;
          }
        }

        if (!found)
        {
          // we didn't find any similar ones
          // so just add it
          this->add_one(junctions, is_reversed);
        }
      }
    }
  }
  else // multi-exon reads
  {
    if (this->junction_dict.size() == 0)
    {
      // there's nothing in the junction dict yet
      // just add it
      this->add_one(junctions, is_reversed);
    }
    else
    {
      // there's stuff in junction_dict
      // check to see if there's an exact match

      if (this->junction_dict.count(junctions.junctions) == 0)
      {
        // the key is not present in the dict yet
        // just add it
        this->add_one(junctions, is_reversed);
      }
      else
      {
        // the key is not directly present in the dict
        // check for similar keys

        bool found = false;

        for (auto & [current_key, current_value] : this->junction_dict)
        {
          if (junctions.junctions.size() == current_key.size())
          {
            // they are the same size. see if they are similar
            
            bool similar = true;
            for (int i = 0; i < current_key.size(); i++)
            {
              if (abs(junctions.junctions[i] - current_key[i]) > this->MAX_DIST)
              {
                // they are too far apart
                similar = false;
                break;
              }
            }

            if (similar)
            {
              // if they are still similar,
              // then we have found an entry that is close enough
              
              this->update_one(junctions, current_key, is_reversed);
              found = true;
              break;
            }
          }
        }
        
        if (!found)
        {
          // we went through and found nothing similar enough
          // so just add it to the dict
          this->add_one(junctions, is_reversed);
        }
      }
    }
  }
};

void Isoforms::add_one(Junctions junctions, bool strand)
{
  if (junctions.junctions.size() == 0) // single-exon
  {
    std::pair<int, int> key = 
    {junctions.left, junctions.right};

    this->single_block_dict[key] = this->single_blocks.size();
    this->single_blocks.push_back({key});
    this->strand_cnt[key] = {};

    int multiplier = 1;
    if (strand) {multiplier = -1;}
    this->strand_cnt[key].push_back(multiplier * this->STRAND_SPECIFIC);
  }
  else // multi-exon reads
  {
    this->junction_dict[junctions.junctions] = this->junction_list.size();
    this->junction_list.push_back({junctions.junctions});
    this->left.push_back(junctions.left);
    this->right.push_back(junctions.right);
    this->lr_pair[junctions.junctions] = {};
    
    int multiplier = 1;
    if (strand) {multiplier = -1;}
    this->lr_pair[junctions.junctions].push_back(multiplier * this->STRAND_SPECIFIC);
  }
}

void Isoforms::update_one(Junctions junctions, std::vector<int> key, bool strand)
{
  if (junctions.junctions.size()==0) // single-exon reads
  {
    std::pair<int, int> new_element = 
    {junctions.left, junctions.right};

    // we need to express the key as a pair here,
    // because we're looking through single_block_dict which is full of pairs
    // honestly this is gonna become a big problem with this 
    std::pair<int, int> key_pair =
    {key[0], key[1]};

    this->single_blocks[this->single_block_dict[key_pair]].push_back(new_element);
  }
  else // multi-exon reads
  {
    this->junction_list[this->junction_dict[key]].push_back(junctions.junctions);
    this->left.push_back(junctions.left);
    this->right.push_back(junctions.right);
    this->lr_pair[key].push_back(std::vector<int> {junctions.left, junctions.right})
  }
}