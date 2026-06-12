// g++ main.cpp -llapacke

#include <lapacke.h>
#include <math.h>
#include <stdio.h>

#define DIM  3
#define DIM2 9

static void F( double io[ DIM ] );
static void jacobian( double x[ DIM ], double o[ DIM2 ] );

static void root_newton_step( double x[ DIM ] )
{
    static double a[ DIM2 ];
    static double b[ DIM ];
    static int ipiv[ DIM ];

    static double fx[ DIM ];
    static double old_fx[ DIM ];
    static double new_x[ DIM ];

    // setup A * X = B  =>  J(x) * (new_x - x) = -F(x)

    jacobian( x, a );
    for ( int i = 0; i < DIM; i++ )
        b[ i ] = x[ i ];
    F( b );
    for ( int i = 0; i < DIM; i++ )
        old_fx[ i ] = b[ i ];
    for ( int i = 0; i < DIM; i++ )
        b[ i ] = -b[ i ];

    int info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, DIM, 1, a, DIM, ipiv, b, 1 );

    // solution stored in b

    if ( info ) {
        printf( "root_newton_step: dgesv failed\n" );
        return;
    }

    double alpha = 1;
backtrack:
    for ( int i = 0; i < DIM; i++ ) {
        new_x[ i ] = x[ i ] + alpha * b[ i ];
        fx[ i ] = new_x[ i ];
    }
    F( fx );

    // keep changes in x small
    // double norm1 = 0;
    // for ( int i = 0; i < DIM; i++ )
    //     norm1 += fmax( norm1, abs( new_x[ i ] - x[ i ] ) );
    // if ( norm1 > 1. ) {
    //     alpha /= 2.;
    //     goto backtrack;
    // }

    // keep changes in fx small
    double norm2 = 0;
    for ( int i = 0; i < DIM; i++ )
        norm2 += fmax( norm2, abs( fx[ i ] - old_fx[ i ] ) );
    if ( norm2 > 1. ) {
        alpha /= 2.;
        printf( "root_newton_step: backtrack\n" );
        goto backtrack;
    }

    for ( int i = 0; i < DIM; i++ )
        x[ i ] = new_x[ i ];
}

static void root_newton( double x[ DIM ] )
{
    static double fx[ DIM ];
    for ( int i = 0; i < 10; i++ ) {
        for ( int k = 0; k < DIM; k++ )
            fx[ k ] = x[ k ];
        F( fx );

        printf( "root_newton: " );
        for ( int k = 0; k < DIM; k++ ) {
            printf( "F(%lf) = %lf ", x[ k ], fx[ k ] );
        }
        printf( "\n" );

        root_newton_step( x );
    }

    printf( "\n" );
}

static double V_A = 10.; // early voltage

static double eber_moll( double vbe, double vce )
{
    return 1E-13 * exp( vbe * 40. ) * ( 1. + vce / V_A );
}

static double eber_moll_partial_vbe( double vbe, double vce )
{
    return 40. * eber_moll( vbe, vce );
}

static double eber_moll_partial_vce( double vbe, double vce )
{
    return 1E-13 * exp( vbe * 40. ) * ( 1. / V_A );
}

static double early_beta( double vce )
{
    double beta = 100;
    return beta * ( 1 + vce / V_A );
}

static double early_beta_partial()
{
    double beta = 100;
    return beta / V_A;
}

static double i1, i2, i3, i4, i5;

static void F( double io[ DIM ] )
{
    double v1 = io[ 0 ];
    double v2 = io[ 1 ];
    double v3 = io[ 2 ];

    i1 = ( 5 - v1 ) / 1000;
    i2 = ( v1 - v3 ) / 1000;
    i3 = ( 9 - v2 ) / 1000;
    i4 = eber_moll( v2 - v3, 0 );
    i5 = eber_moll( v3 - 0, 0 );

    // KCL equations for each node
    io[ 0 ] = i1 - i2;
    io[ 1 ] = i3 - i4;
    io[ 2 ] = i2 + i4 - i5;
}

static void jacobian( double x[ DIM ], double o[ DIM2 ] )
{
    double v1 = x[ 0 ];
    double v2 = x[ 1 ];
    double v3 = x[ 2 ];

    // node 1 KCL derivatives
    o[ 0 ] = -1. / 1000. - ( 1. / 1000. );
    o[ 1 ] = 0. - 0.;
    o[ 2 ] = 0. - ( -1. / 1000. );

    // node 2 KCL derivatives
    o[ 3 ] = 0. - 0.;
    o[ 4 ] = -1. / 1000. - ( eber_moll_partial_vbe( v2 - v3, 0. ) * 1. );
    o[ 5 ] = 0. - ( eber_moll_partial_vbe( v2 - v3, 0. ) * -1. );

    // node 3 KCL derivatives
    o[ 6 ] = ( 1. / 1000. ) + ( 0. ) - ( 0. );
    o[ 7 ] = ( 0. ) + ( eber_moll_partial_vbe( v2 - v3, 0. ) * 1. ) - ( 0. );
    o[ 8 ] = ( -1. / 1000. ) + ( eber_moll_partial_vbe( v2 - v3, 0. ) * -1. ) -
             ( eber_moll_partial_vbe( v3 - 0., 0. ) * 1. );
}

int main()
{
    static double x[ DIM ];
    static double fx[ DIM ];

    x[ 0 ] = 0;
    x[ 1 ] = 0;
    x[ 2 ] = 0;

    root_newton( x );

    for ( int i = 0; i < DIM; i++ )
        printf( "%lf V  ", x[ i ] );
    printf( "\n" );

    printf(
        "%lf mA  %lf mA  %lf mA  %lf mA  %lf mA\n",
        i1 * 1000,
        i2 * 1000,
        i3 * 1000,
        i4 * 1000,
        i5 * 1000
    );

    return 0;
}
