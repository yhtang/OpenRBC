/******************************************************************************
!@ This file is part of the OpenRBC code for protein-resolution red blood cells
!@ Copyright (c) 2016 Yu-Hang Tang, Lu Lu
!@ Licensed under the Apache License, Version 2.0 (the "License")
!@ you may not use this file except in compliance with the License.
!@ You may obtain a copy of the License at
!@             http://www.apache.org/licenses/LICENSE-2.0
!@ Unless required by applicable law or agreed to in writing, software
!@ distributed under the License is distributed on an "AS IS" BASIS,
!@ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!@ See the License for the specific language governing permissions and
!@ limitations under the License.
******************************************************************************/
#ifndef OPENRBC_TRAJECTORY_H_
#define OPENRBC_TRAJECTORY_H_

#include "container.h"
#include "forcefield.h"
#include "voronoi.h"
#include "runtime_parameter.h"
#include <fstream>
#include <cmath>

namespace openrbc {

struct Serializer {
    Serializer( std::ostream & file ) : ost( file ) {}
    ~Serializer() { ost << std::flush; }
    inline void write_title( std::string title ) {
        assert( title.size() <= 8 );
        ost.write( title.c_str(), title.size() );
        for ( int i = 0; i < 8 - title.size(); ++i ) ost.put( 0 );
    }
    template<typename T> friend inline Serializer & operator << ( Serializer & s, T const & var ) {
        static_assert( std::is_pod<T>::value, "Only POD types can be auto-serialized" );
        s.ost.write( ( char const * ) &var, sizeof( T ) );
        return s;
    }
protected:
    std::ostream & ost;
};

struct Deserializer {
    Deserializer( std::istream & file ) : ist( file ) {}
    inline std::string read_title() {
        char buffer[9];
        memset( buffer, 0, 9 );
        for ( int i = 0; i < 8; ++i ) buffer[i] = ist.get();
        return buffer;
    }
    template<typename T> friend inline Deserializer & operator >> ( Deserializer & d, T & var ) {
        static_assert( std::is_pod<T>::value, "Only POD types can be auto-serialized" );
        d.ist.read( ( char * ) &var, sizeof( T ) );
        return d;
    }
protected:
    std::istream & ist;
};

// save the coordinate to file using the ATOM format
template<class CONTAINER1, class CONTAINER2>
void save_frame( std::ostream & file,
                 const CONTAINER1 & model1, const CONTAINER2 & model2,
                 const VCellList & cl_lipid, const VCellList & cl_protein,
                 RTParameter const & param ) {
    Service<Timers>::call()["save_frame"].start();

    Serializer stream( file );

    stream.write_title( "FRAMEBEG" );
    stream << param.nstep ;
    stream.write_title( "NATOM" );
    stream << model1.size() + model2.size();
    stream.write_title( "IDENTITY" );
    for ( std::size_t i = 0 ; i < model1.size() ; ++i ) stream << model1.tag[i] << model1.type[i];
    for ( std::size_t i = 0 ; i < model2.size() ; ++i ) stream << model2.tag[i] << model2.type[i];
    if ( param.dump_field & DumpField::position ) {
        stream.write_title( "POSITION" );
        for ( std::size_t i = 0 ; i < model1.size() ; ++i ) stream << model1.x[i][0] << model1.x[i][1] << model1.x[i][2];
        for ( std::size_t i = 0 ; i < model2.size() ; ++i ) stream << model2.x[i][0] << model2.x[i][1] << model2.x[i][2];
    }
    if ( param.dump_field & DumpField::velocity ) {
        stream.write_title( "VELOCITY" );
        for ( std::size_t i = 0 ; i < model1.size() ; ++i ) stream << model1.v[i][0] << model1.v[i][1] << model1.v[i][2];
        for ( std::size_t i = 0 ; i < model2.size() ; ++i ) stream << model2.v[i][0] << model2.v[i][1] << model2.v[i][2];
    }
    if ( param.dump_field & DumpField::rotation ) {
        stream.write_title( "ROTATION" );
        for ( std::size_t i = 0 ; i < model1.size() ; ++i ) stream << model1.n[i][0] << model1.n[i][1] << model1.n[i][2];
        for ( std::size_t i = 0 ; i < model2.size() ; ++i ) stream << model2.n[i][0] << model2.n[i][1] << model2.n[i][2];
    }
    if ( param.dump_field & DumpField::voronoi ) {
        stream.write_title( "VORONOI" );
        for ( std::size_t i = 0 ; i < model1.size() ; ++i ) stream << cl_lipid.affiliation[i];
        for ( std::size_t i = 0 ; i < model2.size() ; ++i ) stream << cl_protein.affiliation[i];
    }
    if ( param.dump_field & DumpField::force ) {
        stream.write_title( "FORCE" );
        for ( std::size_t i = 0 ; i < model1.size() ; ++i ) stream << model1.f[i][0] << model1.f[i][1] << model1.f[i][2];
        for ( std::size_t i = 0 ; i < model2.size() ; ++i ) stream << model2.f[i][0] << model2.f[i][1] << model2.f[i][2];
    }
    stream.write_title( "FRAMEEND" );

    Service<Timers>::call()["save_frame"].stop();
}

template<class CONTAINER>
bool read_frame( std::istream & file, CONTAINER & model, RTParameter & param ) {
    Service<Timers>::call()["read_frame"].start();

    Deserializer stream( file );

    while ( !file.eof() ) {
        std::string section = stream.read_title();
        if ( file.eof() ) return false;
        if ( section == "FRAMEBEG" ) {
            stream >> param.nstep;
        } else if ( section ==  "FRAMEEND" ) {
            break;
        } else if ( section ==  "NATOM" ) {
            decltype( model.size() ) n;
            stream >> n;
            model.resize( n );
        } else if ( section ==  "IDENTITY" ) {
            for ( std::size_t i = 0 ; i < model.size() ; ++i ) stream >> model.tag[i] >> model.type[i];
        } else if ( section ==  "POSITION" ) {
        	param.dump_field |= DumpField::position;
            for ( std::size_t i = 0 ; i < model.size() ; ++i ) stream >> model.x[i][0] >> model.x[i][1] >> model.x[i][2];
        } else if ( section ==  "VELOCITY" ) {
        	param.dump_field |= DumpField::velocity;
            for ( std::size_t i = 0 ; i < model.size() ; ++i ) stream >> model.v[i][0] >> model.v[i][1] >> model.v[i][2];
        } else if ( section ==  "ROTATION" ) {
        	param.dump_field |= DumpField::rotation;
            for ( std::size_t i = 0 ; i < model.size() ; ++i ) stream >> model.n[i][0] >> model.n[i][1] >> model.n[i][2];
        } else if ( section ==  "VORONOI" ) {
        	param.dump_field |= DumpField::voronoi;
            std::remove_reference<decltype( VCellList::affiliation[0] )>::type dummy;
            for ( std::size_t i = 0 ; i < model.size() ; ++i ) stream >> dummy;
        } else if ( section ==  "FORCE" ) {
        	param.dump_field |= DumpField::force;
            for ( std::size_t i = 0 ; i < model.size() ; ++i ) stream >> model.f[i][0] >> model.f[i][1] >> model.f[i][2];
        } else {
            printf( "Unknown section %s in file\n", section );
        }
    }

    Service<Timers>::call()["read_frame"].stop();

    return true;
}

}

#endif /* TRAJECTORY_H_ */
