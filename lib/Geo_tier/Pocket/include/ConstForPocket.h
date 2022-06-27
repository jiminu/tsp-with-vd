#ifndef _CONST_FOR_POCKET_H
#define _CONST_FOR_POCKET_H


const int NUM_POCKETS = 10;

const float POCKET_DEFAULT_INNER_BETA_VALUE = 1.4f;
const float POCKET_DEFAULT_OUTER_BETA_VALUE = 3.0f;

const float POCKET_COLOR_EXTRANEOUS_ENTITY[3] = {0.7f, 0.7f, 0.7f};
const float POCKET_TRANSPARANCY_EXTRANEOUS_ENTITY = 1.0f;


const float POCKET_COLOR[NUM_POCKETS][3] = {
                    {0.8f,  0.0f,  0.0f},		//����
                    {0.0f,  0.2f,  0.8f},		//�Ķ�
                    {0.0f,  0.4f,  0.0f},		//�ʷ�
                    {1.0f,  1.0f,  0.0f},		//���                  
					{0.8f,  0.2f,  0.6f},		//����
                    {1.0f,  0.4f,  0.0f},		//��Ȳ
                    {0.2f,  0.2f,  0.6f},		//����
                    {0.6f,  0.4f,  0.2f},		//����
                    {0.4f,  0.8f,  1.0f},		//�ϴ�
                    {0.6f,  1.0f,  0.6f},		//����
                };
/*
const float POCKET_COLOR[NUM_POCKETS][3] = {
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},
                    {0.8f,  0.0f,  0.0f},//����                    
                };
*/
const int COMPARE_WRT_NUM_OF_FACES  = 0;
const int COMPARE_WRT_AREA_OF_FACES = 1;

#endif