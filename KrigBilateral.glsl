// KrigBilateral by Shiandow
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library.


//!HOOK CHROMA
//!BIND HOOKED
//!BIND LUMA
//!SAVE LOWRES_Y
//!HEIGHT LUMA.h
//!WHEN CHROMA.w LUMA.w <

// -- Downscaling --
#define offset    (-vec2(0.0, 0.0)*LUMA_size*CHROMA_pt)

#define factor    ((LUMA_pt*CHROMA_size.x)[axis])

#define axis 0

#define Kernel(x) cos(acos(-1.0)*x/2.0)

vec4 hook() {
    // Calculate bounds
    float low  = floor((LUMA_pos - LUMA_pt) * LUMA_size - offset + 0.5)[axis];
    float high = floor((LUMA_pos + LUMA_pt) * LUMA_size - offset + 0.5)[axis];

    float W = 0.0;
    vec4 avg = vec4(0);
    vec2 pos = LUMA_pos;

    for (float k = 0.0; k < high - low; k++) {
        pos[axis] = LUMA_pt[axis] * (k + low + 0.5);
        float rel = (pos[axis] - LUMA_pos[axis])*CHROMA_size[axis] + offset[axis]*factor;
        float w = Kernel(rel);

        vec4 y = textureLod(LUMA_raw, pos, 0.0).xxxx;
        y.y *= y.y;
        avg += w * y;
        W += w;
    }
    avg /= W;

    return avg;
}

//!HOOK CHROMA
//!BIND HOOKED
//!BIND LOWRES_Y
//!SAVE LOWRES_Y
//!WHEN CHROMA.h LUMA.h <

// -- Downscaling --
#define offset    (-vec2(0.0, 0.0)*LOWRES_Y_size*CHROMA_pt)

#define factor    ((LOWRES_Y_pt*CHROMA_size)[axis])

#define axis 1

#define Kernel(x) cos(acos(-1.0)*x/2.0)

vec4 hook() {
    // Calculate bounds
    float low  = floor((LOWRES_Y_pos - LOWRES_Y_pt) * LOWRES_Y_size - offset + 0.5)[axis];
    float high = floor((LOWRES_Y_pos + LOWRES_Y_pt) * LOWRES_Y_size - offset + 0.5)[axis];

    float W = 0.0;
    vec4 avg = vec4(0);
    vec2 pos = LOWRES_Y_pos;

    for (float k = 0.0; k < high - low; k++) {
        pos[axis] = LOWRES_Y_pt[axis] * (k + low + 0.5);
        float rel = (pos[axis] - LOWRES_Y_pos[axis])*CHROMA_size[axis] + offset[axis]*factor;
        float w = Kernel(rel);

        avg += w * textureLod(LOWRES_Y_raw, pos, 0.0);
        W += w;
    }
    avg /= W;

    avg.y = avg.y - avg.x*avg.x;
    return avg;
}

//!HOOK CHROMA
//!BIND HOOKED
//!BIND LUMA
//!BIND LOWRES_Y
//!WIDTH LUMA.w
//!HEIGHT LUMA.h

// -- KrigBilateral --

// -- Convenience --
#define sqr(x)   dot(x,x)
#define bitnoise 1.0/(2.0*255.0)
#define noise    0.05//5.0*bitnoise
#define chromaOffset vec2(0.0, 0.0)

// -- Window Size --
#define taps 2
#define even (float(taps) - 2.0 * floor(float(taps) / 2.0) == 0.0)
#define minX int(1.0-ceil(float(taps)/2.0))
#define maxX int(floor(float(taps)/2.0))

#define Kernel(x) (cos(acos(-1.0)*(x)/float(taps))) // Hann kernel

// -- Input processing --
#define GetY(x,y)  LOWRES_Y_tex(LOWRES_Y_pt*(pos+vec2(x,y)+vec2(0.5))).xy
#define GetUV(x,y) CHROMA_tex(CHROMA_pt*(pos+vec2(x,y)+vec2(0.5))).xy

#define N (taps*taps - 1)

#define M(i,j) Mx[min(i,j)*N + max(i,j) - min(i,j)*(min(i,j)+1)/2]

//#define C(i,j) (1.0 / (1.0 + 0.5*(sqr(X[i].x - X[j].x)/localVar + sqr((coords[i] - coords[j])/radius))) + 0.25 * (X[i].x - y) * (X[j].x - y) / localVar)
//#define c(i)   (1.0 / (1.0 + 0.5*(sqr(X[i].x - y)/localVar + sqr((coords[i] - offset)/radius))))
#define C(i,j) (inversesqrt(1.0 + (X[i].y + X[j].y)/localVar) * exp(-0.5*(sqr(X[i].x - X[j].x)/(localVar + X[i].y + X[j].y) + sqr((coords[i] - coords[j])/radius))) + 0.25 * (X[i].x - y) * (X[j].x - y) / localVar)
#define c(i) (inversesqrt(1.0 + X[i].y/localVar) * exp(-0.5*(sqr(X[i].x - y)/(localVar + X[i].y) + sqr((coords[i] - offset)/radius))))

vec4 hook() {
    float y = LUMA_tex(LUMA_pos).x;

    // Calculate position
    vec2 pos = LUMA_pos * LOWRES_Y_size - chromaOffset - vec2(0.5);
    vec2 offset = pos - (even ? floor(pos) : round(pos));
    pos -= offset;

    vec2 coords[N+1];
    vec4 X[N+1];
    int i=0;
    for (int xx = minX; xx <= maxX; xx++)
    for (int yy = minX; yy <= maxX; yy++)
        if (!(xx == 0 && yy == 0)) {
            coords[i] = vec2(xx, yy);
            X[i++] = vec4(GetY(xx, yy), GetUV(xx, yy));
        }

    coords[N] = vec2(0, 0);
    X[N] = vec4(GetY(0, 0), GetUV(0, 0));

    vec4 total = vec4(0);
    for (int i=0; i<N+1; i++) {
        vec2 w = clamp(1.5 - abs(coords[i] - offset), 0.0, 1.0);
        total += w.x*w.y*vec4(X[i].x, X[i].x*X[i].x, X[i].y, 1.0);
    }
    total.xyz /= total.w;
    float localVar = sqr(noise) + max(0.0, total.y - total.x*total.x + sqr(y - total.x) + total.z);
    float radius = mix(1.0, 2.0, sqr(noise) / localVar);

    float Mx[N*(N+1)/2];
    float b[N];

    for (int i=0; i<N; i++) {
    	b[i] = c(i) - c(N) - ( C(i,N) - C(N,N) ); 
        for (int j=i; j<N; j++) M(i, j)  = C(i,j) - C(j,N) - ( C(i,N) - C(N,N) );
    }

    for (int j=1; j<N; j++) {
        b[j] -= b[0] * M(j, 0) / M(0, 0); 
        for (int k=j; k<N; k++) M(j, k) -= M(0, k) * M(j, 0) / M(0, 0);
    }
    for (int j=2; j<N; j++) {
        b[j] -= b[1] * M(j, 1) / M(1, 1); 
        for (int k=j; k<N; k++) M(j, k) -= M(1, k) * M(j, 1) / M(1, 1);
    }

    vec4 interp = X[N];
    for (int i=N-1; i>=0; i--) {
        for (int j=i+1; j<N; j++) {
            b[i] -= M(i, j) * b[j];
        }
        b[i] /= M(i, i);
        interp += b[i] * (X[i] - X[N]);
    }

    return interp.zwxx;
}
