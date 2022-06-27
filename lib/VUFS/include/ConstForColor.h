#ifndef _CONST_FOR_COLOR_
#define _CONST_FOR_COLOR_

typedef struct {
    double R;
    double G;
    double B;
    double ALPHA;
} Color4DrawingObj;

enum ColorMode { RED, GREEN, DARK_GREEN, BLUE, YELLOW, MAIZE, GRAY, BLACK, WHITE };

const int NUM_OF_COLOR_TYPE = 9;

const Color4DrawingObj DEFAULT_COLOR_OBJ [NUM_OF_COLOR_TYPE] =
{
    { 1.0, 0.0, 0.0, 1.0  },
    { 0.0, 1.0, 0.0, 1.0  },
    { 0.2, 0.8, 0.0, 1.0  },
    { 0.0, 0.0, 1.0, 1.0  },
    { 1.0, 1.0, 0.0, 1.0  },
	{ 0.945, 0.773, 0.282, 1.0  },
    { 0.686275, 0.686275, 0.686275, 1.0  },
    { 0.0, 0.0, 0.0, 1.0  },
    { 1.0, 1.0, 1.0, 1.0  }
};

#endif // _CONST_FOR_COLOR_

