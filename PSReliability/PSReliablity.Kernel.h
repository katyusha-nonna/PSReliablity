#pragma once

#include "PSReliablity.Model.h"
#include "PSReliablity.Config.h"
#include <iomanip>
#include <fstream>
#include <algorithm>

#ifndef  PSRELIABLITY_KERNEL_H
#define PSRELIABLITY_KERNEL_H

// 单个元件生成裕度表
OutageTable* CreateTableForComponent(BaseModel* cp, int start, int length, int duration);

// 并联卷积运算公式
OutageTable ParallelConvolution(OutageTable cp1, OutageTable cp2);
// 串联等效运算公式
OutageTable SeriesEquivalence(OutageTable cp1, OutageTable cp2);

// 发电系统可靠性评估
void GenerationSystemEvaluate(std::ifstream& input, std::ofstream& output);

#endif
