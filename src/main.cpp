#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>
#include <vector>

#include "console_option.hpp"

using namespace std;


using rng_type  = mt19937_64;
using rng_int_t = typename rng_type::result_type;
auto constexpr rng_wsize = rng_type::word_size;


int main( int argc, char** argv )
{
	try
	{
		/**
		***************************************************
		** options.
		***************************************************
		**/
		auto args = nsl::console_option( argc, argv );

		string gen_file = args.read( "-g", "--genotype" );
		string map_file = args.read( "-m", "--map"      );
		string ped_file = args.read( "-p", "--pedigree" ); // gives the number of offsprings.

		size_t seed = random_device{}() | args.read( "-s", "-seed" );


		/**
		***************************************************
		** declaring variables.
		***************************************************
		**/
		// expecting a one column file with distance (in Morgan) between
		// two snps, and a 0 at the first snp of a chromosome.
		auto map          = vector< double >{};

		// expecting a two column file with parent (1-encoded) indices.
		auto pedigree     = vector< size_t >{};

		// expecting a nsnp x 2.nind file with (0/1)s, separated by any
		// one character.
		auto old_genotype = vector< vector< rng_int_t > >{};

		// output genotype separating character.
		auto out_sep = ' ';


		/**
		***************************************************
		** reading files.
		***************************************************
		**/
		// map. -------------------------------------------
		if( auto in = ifstream( map_file ); in.fail() )
			throw runtime_error( "invalid map file : " + map_file );
		else
		{
			auto line = string{};
			while( getline( in, line ) )
				map.push_back( stod(line) );
		}

		// pedigree. --------------------------------------
		if( auto in = ifstream( ped_file ); in.fail() )
			throw runtime_error( "invalid pedigree file : " + ped_file );
		else
		{
			auto line = string{};
			auto sep  = ' ';

			while( getline( in, line ) )
			{
				auto pos = line.find( sep );
				pedigree.push_back( stoi(line.substr(0,pos)) );
				pedigree.push_back( stoi(line.substr(pos+1)) );
			}
		}

		// genotype. --------------------------------------
		if( auto in = ifstream( gen_file ); in.fail() )
			throw runtime_error( "invalid genotype file : " + gen_file );
		else
		{
			auto line = string{};
			auto c64  = size_t{0};

			while( getline( in, line ) )
			{
				old_genotype.resize( (line.size()+1)/2 );

				// the first snp is the highest bit.
				/// @dev : genotypes are 0/1s only, spaced by one character : there are
				/// 		   at even positions.
				if( c64 % rng_wsize != 0 )
					for( auto i = size_t{0}, end = old_genotype.size(); i < end; ++i )
						old_genotype[i].back() = (old_genotype[i].back() << 1) + (line[2*i] == '1');
				else
					for( auto i = size_t{0}, end = old_genotype.size(); i < end; ++i )
						old_genotype[i].push_back( line[2*i] == '1' );

				++c64;
			}
		}

		/**
		***************************************************
		** checking.
		***************************************************
		**/
		if( (old_genotype[0].size() - 1) * rng_wsize + map.size() % rng_wsize != map.size() )
			throw runtime_error( "invalid file format : number of snps does not correspond." );


		/**
		***************************************************
		** assigning variables.
		***************************************************
		**/
		// genotype related variables. --------------------
		auto last_batch   = map.size() % rng_wsize;
		auto nbatch       = old_genotype[0].size();
		auto new_genotype =
			vector< vector < rng_int_t > >( pedigree.size(), vector< rng_int_t >( nbatch ) );

		auto chr_end      = vector< size_t >{};
		for( auto i = size_t{1}, end = map.size(); i < end; ++i )
			if( map[i] == 0 )
				chr_end.push_back(i); // +1 for out of scope.
		chr_end.push_back( map.size() );

		// map related variables. -------------------------
		// transforming %map into a contiguous partial sum.
		auto gap = double{0};
		for( auto i = size_t{1}, end = map.size(); i < end; ++i )
		{
			if( map[i] == 0 )
				gap = map[i-1];

			map[i] += gap;
		}

		auto length = map.back();

		// random variables. ------------------------------
		auto rng = rng_type{ seed };
		auto p   = poisson_distribution< size_t >( length );
		auto u   = uniform_real_distribution< double >( 0, length );


		/**
		***************************************************
		** simulating meioses.
		***************************************************
		**/
		for( auto i = size_t{0}, iend = new_genotype.size(); i < iend; ++i )
		{
			// declaring the mask. ------------------------
			auto mask = vector< rng_int_t >( nbatch );

			auto update_mask = [&]( auto beg, auto end, auto label )
			{
				/// @dev : 0 : first chr; 1 : second chr.
				/// @dev : 0 padding is neutral with xor.
				/// @dev : last batch is accounted for by "shortening" its length.
				/// @dev : adapted for [first;last).

				// infering batch from snp position.
				/// @dev : %end is out of the range.
				auto batch1 = beg / rng_wsize;
				auto batch2 = end / rng_wsize;
				auto snp1   = beg % rng_wsize;
				auto snp2   = end % rng_wsize;

				// updating the mask.
				mask[batch1] ^= (-label) >> snp1 + ( batch1 != nbatch - 1 ? 0 : rng_wsize - last_batch );
				for( auto s = batch1 + 1; s < batch2; ++s )
					mask[s]  ^= -label;
				mask[batch2] ^= (-label) << ( batch2 != nbatch - 1 ? rng_wsize : last_batch ) - snp2;
			};

			// crossovers ---------------------------------
			// generating crossovers (Haldane's model of recombination).
			auto crossovers = vector< double >( p(rng) );
			for( auto& cross : crossovers ) cross = u(rng);
			sort( crossovers.begin(), crossovers.end() );

			// constructing the mask for recombinations.
			{
				auto it  = map.begin();
				auto beg = size_t{0};

				for( auto c = size_t{0}, cend = crossovers.size(); c < cend; ++c )
				{
					// infering snp positions.
					it       = lower_bound( it, map.end(), crossovers[c] );
					auto end = distance( map.begin(), it );

					// generating the label.
					auto label = rng_int_t{ c & 1 };

					// updating the mask.
					update_mask( beg, end, label );

					// updating the range iterators.
					beg = end;
				}
			}

			// mendelian sampling. ------------------------
			// constructing the mask for mendelian sampling.
			{
				auto beg = size_t{0};

				for( auto end : chr_end )
				{
					// generating the label.
					/// @dev : 0 is no change, 1 is change.
					auto label = rng_int_t{ rng() & 1 };

					// updating the mask.
					update_mask( beg, end, label );

					// updating the range iterators.
					beg = end;
				}
			}

			// masking. -----------------------------------
			auto& new_gamete  = new_genotype[i];
			auto& old_chr_1   = old_genotype[pedigree[i]*2-2];
			auto& old_chr_2   = old_genotype[pedigree[i]*2-1];
			for( auto s = size_t{0}; s < nbatch; ++s )
				new_gamete[s] = ( old_chr_1[s] & ~mask[s] ) | ( old_chr_2[s] & mask[s] );
		}


		/**
		***************************************************
		** outputting.
		***************************************************
		**/
		for( auto s = size_t{0}, send = nbatch; s < send; ++s )
		{
			auto byte_size = s != nbatch - 1 ? rng_wsize : last_batch;

			for( auto b = size_t{0}; b < byte_size; ++b )
				for( auto i = size_t{0}, iend = new_genotype.size(); i < iend; ++i )
					cout << ( new_genotype[i][s] >> byte_size - b - 1 & 1 ) << ( i != iend - 1 ? out_sep : '\n' );
		}
	}
	catch( exception& e )
	{
		cerr << e.what() << endl;
		return 1;
	}

	return 0;
}
