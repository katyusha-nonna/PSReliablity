#include "PSReliablity.Kernel.h"

OutageTable* CreateTableForComponent(BaseModel* cp, int start, int length, int duration)
{
	// ��Ž��
	auto result = new OutageTable(start);
	result->exactProbability.resize(length);
	result->incrementalFrequency.resize(length);
	result->accumulatedProbability.resize(length);
	result->accumulatedFrequency.resize(length);
	// �����
	// �������ǿ��ͣ����q�����и���p
	if (auto gen = dynamic_cast<GeneratorModel*>(cp))
	{
		// ���ڷ����Ԫ��
		// ���ȼ���ǿ��ͣ���ʺ����и���
		double q = gen->lambda / (gen->lambda + gen->miu);
		double p = 1 - q;
		// Ȼ���������Ƶ��
		double f = q * gen->miu;
		// Ȼ������ͣ�˱�
		result->exactProbability[0] = p;
		result->exactProbability[length - 1] = q;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedProbability[i] = q;
		}
		result->accumulatedProbability[0] = 1;
		result->incrementalFrequency[0] = -f;
		result->incrementalFrequency[length - 1] = f;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedFrequency[i] = f;
		}
		result->accumulatedFrequency[0] = 0;
	}
	else if (auto load = dynamic_cast<LoadModel*>(cp))
	{
		// ���ڸ���Ԫ��
		// ͳ�Ƹ�����ֵ���ִ���������ȷ�и��ʺ��ۼƸ���
		for(int i = 0; i < length; i++)
		{
			double num = std::count(load->loadCurve.begin(), load->loadCurve.end(), i);
			result->exactProbability[i] = num / duration;
		}
		result->accumulatedProbability[length-1] = result->exactProbability[length - 1];
		for(int i = length - 2; i > -1; i--)
		{
			result->accumulatedProbability[i] = result->accumulatedProbability[i + 1] + result->exactProbability[i];
		}
		// ɨ�踺�����ߣ���������Ƶ�ʺ��ۼ�Ƶ��
		auto update = [&](int pre, int cur) {
			if (pre < cur)
			{
				// ˵����������
				for (int i = pre; i < cur; i++)
				{
					result->incrementalFrequency[i] += -1;
				}
			}
			else if (pre > cur)
			{
				// ˵�����ɼ���
				for (int i = pre; i > cur; i--)
				{
					result->incrementalFrequency[i] += 1;
				}
			}
			else
			{
				// ����δ�����仯����ʱ�����в���
			}
		};
		// �ȴ����һ��ʱ�̵ĸ���
		update(*(load->loadCurve.rbegin()), *(load->loadCurve.begin()));
		// Ȼ����ʣ�µ�
		for (int i = 1; i < load->loadCurve.size(); i++)
		{
			update(load->loadCurve[i - 1], load->loadCurve[i]);
		}
		// �����ۻ�Ƶ��
		result->accumulatedFrequency[length - 1] = result->incrementalFrequency[length - 1];
		for (int i = length - 2; i > -1; i--)
		{
			result->accumulatedFrequency[i] = result->accumulatedFrequency[i + 1] + result->incrementalFrequency[i];
		}
	}
	else if (auto line = dynamic_cast<LineModel*>(cp))
	{
		// ������·
		// ���ȼ���ǿ��ͣ���ʺ����и���
		double q = line->lambda / (line->lambda + line->miu);
		double p = 1 - q;
		// Ȼ���������Ƶ��
		double f = q * line->miu;
		// Ȼ������ͣ�˱�
		result->exactProbability[0] = p;
		result->exactProbability[length - 1] = q;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedProbability[i] = q;
		}
		result->accumulatedProbability[0] = 1;
		result->incrementalFrequency[0] = -f;
		result->incrementalFrequency[length - 1] = f;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedFrequency[i] = f;
		}
		result->accumulatedFrequency[0] = 0;
	}
	else if (auto tran = dynamic_cast<TransformerModel*>(cp))
	{
		// ���ڱ�ѹ��
				// ���ȼ���ǿ��ͣ���ʺ����и���
				double q = tran->lambda / (tran->lambda + tran->miu);
				double p = 1 - q;
				// Ȼ���������Ƶ��
				double f = q * tran->miu;
				// Ȼ������ͣ�˱�
				result->exactProbability[0] = p;
				result->exactProbability[length - 1] = q;
				for (int i = 1; i < length; i++)
				{
					result->accumulatedProbability[i] = q;
				}
				result->accumulatedProbability[0] = 1;
				result->incrementalFrequency[0] = -f;
				result->incrementalFrequency[length - 1] = f;
				for (int i = 1; i < length; i++)
				{
					result->accumulatedFrequency[i] = f;
				}
				result->accumulatedFrequency[0] = 0;
	}
	return result;
}

OutageTable ParallelConvolution(OutageTable cp1, OutageTable cp2)
{
	// ���ս����ά�ȡ���ʼʱ��
	int length = 0, start = 0;
	length = cp1.accumulatedProbability.size() + cp2.accumulatedProbability.size() -1;
	start = cp1.start + cp2.start;
	// ��Ž��
	OutageTable result(start);
	std::vector<double> exactProbability(length);
	std::vector<double> accumulatedProbability(length);
	std::vector<double> incrementalFrequency(length);
	std::vector<double> accumulatedFrequency(length);
	// lambda������;�����о������
	auto convolution = [&](const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& r) {
		int xSize = x.size(), ySize = y.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < xSize && j <= i; j++)
			{
				if (i - j < ySize) r[i] += x[j] * y[i - j];
			}
		}			
	};
	// lambda������;���Ľ��������(���ȷ�и��ʲ���ľ������λ��y��λ��)
	auto convolution2 = [&](const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& r) {
		int xSize = x.size(), ySize = y.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < xSize; j++)
			{
				if (i - j < ySize)
				{
					if (i - j < 0) r[i] += x[j];
					else r[i] += x[j] * y[i - j];
				}
			}
		}
	};
	// ����ȷ�и��� = cp1��ȷ�и��ʺ�cp2��ȷ�и��ʵľ��
	convolution(cp1.exactProbability, cp2.exactProbability, exactProbability);
	// �����ۻ����� = cp1��ȷ�и��ʺ�cp2���ۻ����ʵľ��
	convolution2(cp1.exactProbability, cp2.accumulatedProbability, accumulatedProbability);
	// ��������Ƶ�� = cp1��ȷ�и��ʺ�cp2������Ƶ�ʵľ�� + cp1������Ƶ�ʺ�cp2��ȷ�и��ʵľ��
	convolution(cp1.exactProbability, cp2.incrementalFrequency, incrementalFrequency);
	convolution(cp1.incrementalFrequency, cp2.exactProbability, incrementalFrequency);
	// �����ۼ�Ƶ�� = cp1��ȷ�и��ʺ�cp2���ۼ�Ƶ�ʵľ�� + cp1������Ƶ�ʺ�cp2���ۼƸ��ʵľ��
	convolution(cp1.exactProbability, cp2.accumulatedFrequency, accumulatedFrequency);
	convolution2(cp1.incrementalFrequency, cp2.accumulatedProbability, accumulatedFrequency);
	// ���ؽ��
	result.exactProbability = std::move(exactProbability);
	result.accumulatedProbability = std::move(accumulatedProbability);
	result.incrementalFrequency = std::move(incrementalFrequency);
	result.accumulatedFrequency = std::move(accumulatedFrequency);
	return std::move(result);
}

OutageTable SeriesEquivalence(OutageTable cp1, OutageTable cp2)
{
	// ���ս����ά�ȡ���ʼʱ��
	int length = 0, start = 0;
	// �����ж�cp1��cp2�Ƿ�ά����ͬ(������ͬ)
	if (cp1.start == cp2.start)
	{
		// ������ͬ��ֱ�Ӽ���
		length = cp1.exactProbability.size();
	}
	else
	{
		// ������ͬ������ֱ�Ӽ��㣬Ҫ���и���
		if (cp1.start < cp2.start)
		{
			length = cp1.exactProbability.size();
			int index = cp2.start - cp1.start;
			std::vector<double> temp;
			temp.assign(cp2.accumulatedProbability.begin() + index, cp2.accumulatedProbability.end());
			cp2.accumulatedProbability = std::move(temp);
			cp2.accumulatedProbability[0] = 1;
			temp.assign(cp2.accumulatedFrequency.begin() + index, cp2.accumulatedFrequency.end());
			cp2.accumulatedFrequency = std::move(temp);
			cp2.accumulatedFrequency[0] = 0;
		}
		else
		{
			length = cp2.exactProbability.size();
			int index = cp1.start - cp2.start;
			std::vector<double> temp;
			temp.assign(cp1.accumulatedProbability.begin() + index, cp1.accumulatedProbability.end());
			cp1.accumulatedProbability = std::move(temp);
			cp1.accumulatedProbability[0] = 1;
			temp.assign(cp1.accumulatedFrequency.begin() + index, cp1.accumulatedFrequency.end());
			cp1.accumulatedFrequency = std::move(temp);
			cp1.accumulatedFrequency[0] = 0;
		}
	}
	start = std::min(cp1.start, cp2.start);
	// ��Ž��
	OutageTable result(start);
	std::vector<double> exactProbability(length);
	std::vector<double> accumulatedProbability(length);
	std::vector<double> incrementalFrequency(length);
	std::vector<double> accumulatedFrequency(length);
	// �����ۼƸ���
	for (int i = 0; i < length; i++)
	{
		accumulatedProbability[i] += cp1.accumulatedProbability[i] + cp2.accumulatedProbability[i] - cp1.accumulatedProbability[i] * cp2.accumulatedProbability[i];
	}
	// �����ۼ�Ƶ��
	for (int i = 0; i < length; i++)
	{
		accumulatedFrequency[i] += cp1.accumulatedFrequency[i] * (1 - cp2.accumulatedProbability[i]);
		accumulatedFrequency[i] += cp2.accumulatedFrequency[i] * (1 - cp1.accumulatedProbability[i]);
	}
	// ����ȷ�и��ʺ�����Ƶ��
	exactProbability[length - 1] = accumulatedProbability[length - 1];
	for (int i = length - 2; i > -1; i--)
	{
		exactProbability[i] = accumulatedProbability[i] - accumulatedProbability[i + 1];
	}
	incrementalFrequency[length - 1] = accumulatedFrequency[length - 1];
	for (int i = length - 2; i > -1; i--)
	{
		incrementalFrequency[i] = accumulatedFrequency[i] - accumulatedFrequency[i + 1];
	}
	// ���ؽ��
	result.exactProbability = std::move(exactProbability);
	result.accumulatedProbability = std::move(accumulatedProbability);
	result.incrementalFrequency = std::move(incrementalFrequency);
	result.accumulatedFrequency = std::move(accumulatedFrequency);
	return std::move(result);
}

void GenerationSystemEvaluate(std::ifstream& input, std::ofstream& output)
{
	// �������ڴ洢Ԫ����Ϣ�ͼ����м����ı���
	// ���ʽ�������ڴ�ż���·��
	ExprTree<OutageTable> systemGraph;
	// �ֵ䣬���ڴ��Ԫ���б�
	std::map<std::string, BaseModel* >  components;
	// ��������
	ReliablityEvaluateConfig config;

	// ��Ҫ�ĳ�ʼ��
	{
		auto comp = [](std::string a, std::string b)->bool {
			if (a == "+")
			{
				if (b == "*")
				{
					return false;
				}
				else
				{
					return true;
				}
			}
			else if (a == "*")
			{
				return true;
			}
			else
			{
				return false;
			}
		};
		systemGraph.AddComp(comp);
		systemGraph.AddOperator("+", ParallelConvolution);
		systemGraph.AddOperator("*", SeriesEquivalence);
	}

	// ��ȡģʽ�Ͳ����Լ�(���ɵ�)ʱ��߶�
	input >> config.evaluateInfo.mode >> config.evaluateInfo.step >> config.evaluateInfo.duration;
	// ��ȡϵͳͳ����Ϣ
	input >> config.systemInfo.numComp;
	input >> config.systemInfo.numGen;
	input >> config.systemInfo.numLoad;
	input >> config.systemInfo.numLine;
	input >> config.systemInfo.numTran;
	// ��ȡ�����
	for (int i = 0; i < config.systemInfo.numGen; i++)
	{
		auto gen = new GeneratorModel();
		input >> gen->name >> gen->ID >> gen->capacity >> gen->lambda >> gen->miu;
		components.try_emplace(gen->ID, gen);
	}
	// ��ȡ����
	for (int i = 0; i < config.systemInfo.numLoad; i++)
	{
		auto load = new LoadModel();
		input >> load->name >> load->ID >> load->capacity;
		load->loadCurve.resize(config.evaluateInfo.duration);
		for (int i = 0; i < config.evaluateInfo.duration; i++)
		{
			input >> load->loadCurve[i];
			load->loadCurve[i] /= config.evaluateInfo.step;
		}
		components.try_emplace(load->ID, load);
	}
	// ��ȡ��·
	for (int i = 0; i < config.systemInfo.numLine; i++)
	{
		auto line = new LineModel();
		input >> line->name >> line->ID >> line->capacity >> line->lambda >> line->miu;
		components.try_emplace(line->ID, line);
	}
	// ��ȡ��ѹ��
	for (int i = 0; i < config.systemInfo.numTran; i++)
	{
		auto tran = new TransformerModel();
		input >> tran->name >> tran->ID >> tran->capacity >> tran->lambda >> tran->miu;
		components.try_emplace(tran->ID, tran);
	}
	// Ϊ����Ԫ���γ�ԣ�ȱ�����
	{
		for (auto& cp : components)
		{
			if (auto gen = dynamic_cast<GeneratorModel*>(cp.second))
			{
				int start = gen->capacity / config.evaluateInfo.step;
				int length = start + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
			else if (auto load = dynamic_cast<LoadModel*>(cp.second))
			{
				int start = 0;
				int length = load->capacity / config.evaluateInfo.step + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
			else if (auto line = dynamic_cast<LineModel* >(cp.second))
			{
				int start = line->capacity / config.evaluateInfo.step;
				int length = start + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
			else if (auto tran = dynamic_cast<TransformerModel*>(cp.second))
			{
				int start = tran->capacity / config.evaluateInfo.step;
				int length = start + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
		}
	}
	// ���ٱ���
	for (auto& cp : components)
	{
		delete cp.second;
		cp.second = nullptr;
	}
	// ��ȡ�������˲��γɱ��ʽ��
	{
		std::string s;
		std::stringstream ss;
		while (s == "")
		{
			getline(input, s);
		}
		ss << s;
		systemGraph.BuildFromInfixExpr(ss);
	}
	auto result = systemGraph.Calcute();
	// ������
	{
		// �������ϵͳͣ�˱�/ԣ�ȱ�
		int TableLength = 25;
		int count = 0;
		int length = result.exactProbability.size();
		output << "���" << std::setw(TableLength) << "ԣ��" << std::setw(TableLength) << "P*" << std::setw(TableLength) << "F*" << std::setw(TableLength)
			<< "P" << std::endl;
		for (int count = 0; count < length; count++)
		{
			output << count << std::setw(TableLength) << (result.start-count) * config.evaluateInfo.step << std::setw(TableLength) << result.accumulatedProbability[count] << std::setw(TableLength) << result.accumulatedFrequency[count] << std::setw(TableLength)
				<< result.exactProbability[count] << std::endl;
		}
		// ����ɿ���ָ��
		double lolp = 0, eens = 0, lolf = 0;
		int k = result.start + 1;
		lolp = k < result.accumulatedProbability.size() ? result.accumulatedProbability[k] : 0;
		for (int i = k; i < length; i++)
		{
			eens += result.accumulatedProbability[i];
		}
		eens *= config.evaluateInfo.step*config.evaluateInfo.duration;
		lolf = k < result.accumulatedFrequency.size() ? result.accumulatedFrequency[k] : 0;
		output << std::endl;
		output << "LOLP" << std::setw(TableLength) << lolp << std::endl;
		output << "EENS" << std::setw(TableLength) << eens << std::endl;
		output << "LOLF" << std::setw(TableLength) << lolf << std::endl;
	}
}