#pragma once

#include "PSReliablity.Model.h"
#include "PSReliablity.Config.h"
#include <iomanip>
#include <fstream>
#include <algorithm>

#ifndef  PSRELIABLITY_KERNEL_H
#define PSRELIABLITY_KERNEL_H

// ����Ԫ������ԣ�ȱ�
OutageTable* CreateTableForComponent(BaseModel* cp, int start, int length, int duration);

// ����������㹫ʽ
OutageTable ParallelConvolution(OutageTable cp1, OutageTable cp2);
// ������Ч���㹫ʽ
OutageTable SeriesEquivalence(OutageTable cp1, OutageTable cp2);

// ����ϵͳ�ɿ�������
void GenerationSystemEvaluate(std::ifstream& input, std::ofstream& output);

#endif
