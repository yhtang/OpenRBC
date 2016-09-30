#include "trajectory.h"
#include <map>
#include <string>
#include <tuple>
#include <fstream>
#include <algorithm>

static std::map<std::string, int ( * )( int, char ** )> commands;

using namespace openrbc;
using namespace openrbc::config;

int main( int argc, char ** argv ) {
    if ( argc < 2 || commands.find( argv[1] ) == commands.end() ) {
        return commands[ "help" ]( argc, argv );
    }
    return commands[ argv[1] ]( argc, argv );
}

struct helper {
    static int action( int argc, char ** argv ) {
        printf( "List of commands:\n" );
        for ( auto & c : commands ) printf( "%s\n", c.first.c_str() );
        return 0;
    }
};
static auto register_helper = commands.emplace( "help", helper::action );

/*-----------------------------------------------------------------------------
Convert ORBC trajectory to LAMMPSTRJ
-----------------------------------------------------------------------------*/
struct convert {
    static int action( int argc, char ** argv ) {
        if ( argc < 4 ) {
            printf( "Usage: %s %s input_filename output_filename\n", argv[0], argv[1] );
            return -1;
        }
        auto from = get_suffix( argv[2] );
        auto to   = get_suffix( argv[3] );
        if ( from == "orbc" && to == "lammpstrj" ) orbc2lammpstrj( argv[2], argv[3] );
        if ( from == "orbc" && to == "namdbin" ) orbc2namdbin( argv[2], argv[3] );

        return 0;
    }
protected:
    static std::string get_suffix( std::string name ) {
        std::string suffix;
        for ( int i = name.size() - 1; i >= 0 && name[i] != '.'; --i ) suffix.insert( 0, 1, name[i] );
        return suffix;
    }

    static int orbc2lammpstrj( std::string from, std::string to ) {
        using namespace config;

        std::ifstream fin( from );
        std::ofstream fout( to );

        RTParameter param( 0, NULL );
        ProteContainer cell( "cell" );
        decltype( VCellList::affiliation ) affiliation;
        int field = DumpField::position + DumpField::voronoi;

        while ( !fin.eof() ) {

            // read ORBC frame

            bool loaded = true;

            while ( !fin.eof() ) {

                Deserializer stream( fin );
                std::string section = stream.read_title();
                if ( fin.eof() ) {
                    loaded = false;
                    break;
                }
                if ( section == "FRAMEBEG" ) {
                    stream >> param.nstep;
                    printf( "Loading frame %d\n", param.nstep );
                } else if ( section ==  "FRAMEEND" ) {
                    loaded = true;
                    break;
                } else if ( section ==  "NATOM" ) {
                    decltype( cell.size() ) n;
                    stream >> n;
                    cell.resize( n );
                    affiliation.resize( n );
                } else if ( section ==  "IDENTITY" ) {
                    for ( std::size_t i = 0 ; i < cell.size() ; ++i ) stream >> cell.tag[i] >> cell.type[i];
                } else if ( section ==  "POSITION" ) {
                    for ( std::size_t i = 0 ; i < cell.size() ; ++i ) stream >> cell.x[i][0] >> cell.x[i][1] >> cell.x[i][2];
                } else if ( section ==  "VELOCITY" ) {
                    for ( std::size_t i = 0 ; i < cell.size() ; ++i ) stream >> cell.v[i][0] >> cell.v[i][1] >> cell.v[i][2];
                } else if ( section ==  "ROTATION" ) {
                    for ( std::size_t i = 0 ; i < cell.size() ; ++i ) stream >> cell.n[i][0] >> cell.n[i][1] >> cell.n[i][2];
                } else if ( section ==  "VORONOI" ) {
                    for ( std::size_t i = 0 ; i < cell.size() ; ++i ) stream >> affiliation[i];
                } else if ( section ==  "FORCE" ) {
                    for ( std::size_t i = 0 ; i < cell.size() ; ++i ) stream >> cell.f[i][0] >> cell.f[i][1] >> cell.f[i][2];
                } else {
                    printf( "Unknown section %s in file\n", section );
                }
            }

            // write LAMMPSTRJ frame

            if ( loaded ) {
                std::vector<int> key( cell.size() );
                for ( int i = 0; i < key.size(); ++i ) key[i] = i;
                std::sort( key.begin(), key.end(), [&cell]( int i, int j ) { return cell.tag[i] < cell.tag[j]; } );

                printf( "Writing frame %d\n", param.nstep );

                // save the coordinate to fout using the ATOM format
                fout << "ITEM: TIMESTEP\n";
                fout << param.nstep << '\n';

                // Number of atoms
                fout << "ITEM: NUMBER OF ATOMS\n";
                fout << cell.size() << '\n';

                // Box Bounds
                fout << "ITEM: BOX BOUNDS ff ff ff\n";
                fout << param.box[0][0] << ' ' << param.box[0][1] << '\n'; // [xmin, xmax]
                fout << param.box[1][0] << ' ' << param.box[1][1] << '\n'; // [ymix, ymax]
                fout << param.box[2][0] << ' ' << param.box[2][1] << '\n'; // [zmin, zmax]

                // output AtrajectoryOMS
                fout << "ITEM: ATOMS id type";
                if ( field & DumpField::position ) fout << " xu yu zu";
                if ( field & DumpField::velocity ) fout << " vx vy vz";
                if ( field & DumpField::rotation ) fout << " mux muy muz";
                if ( field & DumpField::voronoi  ) fout << " vx";
                if ( field & DumpField::force    ) fout << " fx fy fz";
                fout << '\n';

                for ( std::size_t itr = 0 ; itr < cell.size() ; ++itr ) {
                    // (atom's nbr), (atm type), (x,y,z) pos, (x, y, z) mu, (voronoi cell index), (x,y,z) vel, (x,y,z) force of the (i^th) atm
                    int i = key[itr];
                    fout << cell.tag[i] << ' ' << cell.type[i] << ' ';
                    if ( field & DumpField::position ) fout << cell.x[i][0] << ' ' << cell.x[i][1] << ' ' << cell.x[i][2] << ' ';
                    if ( field & DumpField::velocity ) fout << cell.v[i][0] << ' ' << cell.v[i][1] << ' ' << cell.v[i][2] << ' ';
                    if ( field & DumpField::rotation ) fout << cell.n[i][0] << ' ' << cell.n[i][1] << ' ' << cell.n[i][2] << ' ';
                    if ( field & DumpField::voronoi  ) fout << affiliation[i] << ' ';
                    if ( field & DumpField::force    ) fout << cell.f[i][0] << ' ' << cell.f[i][1] << ' ' << cell.f[i][2] << ' ';
                    fout << '\n';
                }
                fout << std::flush;
            }
        }

        return 0;
    }

    static int orbc2namdbin( std::string from, std::string to ) {
        using namespace config;

        std::ifstream fin( from );

        RTParameter param( 0, NULL );
        ProteContainer cell( "cell" );
        decltype( VCellList::affiliation ) affiliation;
        int field = DumpField::position;
        decltype( cell.size() ) effective_size; // due to cleanup of stray lipid

        while ( !fin.eof() ) {

            // read ORBC frame

            bool loaded = true;

            while ( !fin.eof() ) {

                Deserializer stream( fin );
                std::string section = stream.read_title();
                if ( fin.eof() ) {
                    loaded = false;
                    break;
                }
                if ( section == "FRAMEBEG" ) {
                    stream >> param.nstep;
                    printf( "Loading frame %d\n", param.nstep );
                } else if ( section ==  "FRAMEEND" ) {
                    loaded = true;
                    break;
                } else if ( section ==  "NATOM" ) {
                    stream >> effective_size;
                    if ( cell.size() < effective_size ) cell.resize( effective_size );
                    if ( affiliation.size() < effective_size ) affiliation.resize( effective_size );
                } else if ( section ==  "IDENTITY" ) {
                    for ( std::size_t i = 0 ; i < effective_size ; ++i ) stream >> cell.tag[i] >> cell.type[i];
                    for ( std::size_t i = effective_size ; i < cell.size() ; ++i ) cell.tag[i] = cell.size() + 1;
                } else if ( section ==  "POSITION" ) {
                    for ( std::size_t i = 0 ; i < effective_size ; ++i ) stream >> cell.x[i][0] >> cell.x[i][1] >> cell.x[i][2];
                    for ( std::size_t i = effective_size ; i < cell.size() ; ++i ) cell.x[i] = 0;
                } else if ( section ==  "VELOCITY" ) {
                    for ( std::size_t i = 0 ; i < effective_size ; ++i ) stream >> cell.v[i][0] >> cell.v[i][1] >> cell.v[i][2];
                } else if ( section ==  "ROTATION" ) {
                    for ( std::size_t i = 0 ; i < effective_size ; ++i ) stream >> cell.n[i][0] >> cell.n[i][1] >> cell.n[i][2];
                } else if ( section ==  "VORONOI" ) {
                    for ( std::size_t i = 0 ; i < effective_size ; ++i ) stream >> affiliation[i];
                } else if ( section ==  "FORCE" ) {
                    for ( std::size_t i = 0 ; i < effective_size ; ++i ) stream >> cell.f[i][0] >> cell.f[i][1] >> cell.f[i][2];
                } else {
                    printf( "Unknown section %s in file\n", section );
                }
            }

            // write NAMDBIN frame

            if ( loaded ) {
                char filename[256];
                snprintf( filename, 256, to.c_str(), param.nstep );

                std::ofstream fout( filename );
                Serializer stream( fout );

                std::vector<int> key( cell.size() );
                for ( int i = 0; i < key.size(); ++i ) key[i] = i;
                std::sort( key.begin(), key.end(), [&cell]( int i, int j ) { return cell.tag[i] < cell.tag[j]; } );

                printf( "Writing frame %d\n", param.nstep );
                std::uint32_t n = cell.size();
                stream << n;
                for ( std::size_t itr = 0 ; itr < cell.size() ; ++itr ) {
                    int i = key[itr];
                    double x = cell.x[i][0], y = cell.x[i][1], z = cell.x[i][2];
                    stream << x << y << z;
                }
            }
        }

        return 0;
    }
};
static auto register_convert = commands.emplace( "convert", convert::action );
