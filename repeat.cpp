#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
using namespace std;

// �����ַ����Ļ�����ת����
string get_reverse_complement(const string& s) {
    string res = s;
    reverse(res.begin(), res.end());
    for (char& c : res) {
        if (c == 'A') c = 'T';
        else if (c == 'T') c = 'A';
        else if (c == 'C') c = 'G';
        else if (c == 'G') c = 'C';
    }
    return res;
}

int main() {
    string reference, query;
    cout << "Please input the reference sequence:" << endl;
    cin >> reference;
    cout << "Please input the query sequence:" << endl;
    cin >> query;

    int n = query.length();
    int m = reference.length();

    // Ԥ���� reference �������Ӵ����以����ת����
    unordered_map<string, pair<int, bool>> substring_map; // �Ӵ� -> {��ʼλ��, �Ƿ�ת}
    for (int i = 0; i < m; i++) {
        string sub = "";
        for (int j = i; j < m; j++) {
            sub += reference[j];
            substring_map[sub] = {i, false}; // ���Ϊ�ǻ�����ת
        }
    }
    
    // ���� reference �Ļ�����ת�汾
    string reference_rc = get_reverse_complement(reference);
    
    // ���� reference_rc �������Ӵ�
    for (int i = 0; i < m; i++) {
        string sub = "";
        for (int j = i; j < m; j++) {
            sub += reference_rc[j];
            int len = j - i + 1;
            int original_start = m - len - i; // ת��Ϊ reference �е���ʼ�±�
            if (substring_map.find(sub) == substring_map.end()) {
                substring_map[sub] = {original_start, true}; // ���Ϊ������ת
            }
        }
    }

    // ��̬�滮�ָ� query
    vector<int> dp(n + 1, n + 1); // ��ʼ��Ϊ���ֵ
    dp[n] = 0; // �մ�Ϊ 0
    // trace[i] ��¼��λ�� i ����ƥ�����һ��λ�õ���Ϣ��{next_pos, {start_pos, is_inversion}}
    vector<pair<int, pair<int, bool>>> trace(n + 1, {-1, {-1, false}});

    for (int i = n - 1; i >= 0; i--) {
        string sub = "";
        for (int j = i; j < n; j++) {
            sub += query[j];
            if (substring_map.find(sub) != substring_map.end()) {
                int next_pos = j + 1;
                if (dp[i] > 1 + dp[next_pos]) {
                    dp[i] = 1 + dp[next_pos];
                    int start = substring_map[sub].first;
                    bool is_rc = substring_map[sub].second;
                    trace[i] = {next_pos, {start, is_rc}};
                }
            }
        }
    }

    // ���� trace ����ÿ���ֶεĽ����Ϣ���洢Ϊ {start_pos, length, is_inversion}
    // ���� length Ϊÿ���ֶΣ��ظ���λ���ĳ���
    vector<pair<pair<int,int>, bool>> segments;
    int pos = 0;
    while (pos < n) {
        int next = trace[pos].first;
        int start = trace[pos].second.first;
        bool is_rc = trace[pos].second.second;
        int len = next - pos;  // �ظ���λ�ĳ���
        segments.push_back({{start, start + len - 1}, is_rc});
        pos = next;
    }

    // ������������ͬ����� segments ���з���ϲ�
    // �ж�����Ϊ����ʼλ�úͷ�ת��־��ͬ���ظ�������Ϊ������ͬ�ֶεĸ�����
    // ��ÿ���ֶε� length �ǵ�λ����
    struct Output {
        int start_pos;
        int length;
        int repeat_count;
        bool is_inversion;
    };
    vector<Output> outputs;
    for (size_t i = 0; i < segments.size();) {
        int current_start = segments[i].first.first;
        // ��λ���ȼ���Ϊ�������Ҷ˵� - ��ʼ + 1��
        int unit_length = segments[i].first.second - segments[i].first.first + 1;
        bool current_rc = segments[i].second;
        int count = 1;
        size_t j = i + 1;
        while (j < segments.size() && segments[j].first.first == current_start && segments[j].second == current_rc) {
            count++;
            j++;
        }
        outputs.push_back({current_start, unit_length, count, current_rc});
        i = j;
    }

    cout << "Notice:" << endl;
    cout << "1. �������ʽ�У���ǰ�����˳�򿼲�query�����е�ÿһ�ηֱ�������reference�ĺ�λ�ã����ƶ��ٴΣ��Ƿ�ת����Ϣ��"
    cout << "2. ÿ����¼�У�start_pos ��ʾ query �ĸö������� reference �����ʼλ�ã�0 Ϊ��ʼ����λ�ã���length Ϊ�ظ���λ���ȣ�repeat_count Ϊ�ظ����������� 1 �Σ����䱾����is_inversion ��ʾ�ö������Ƿ�ת��" << endl;
    
    // ����������ʽΪ��
    // {'start_pos': X, 'length': Y, 'repeat_count': Z, 'is_inversion': True/False}
    for (auto& out : outputs) {
        cout << "{'start_pos': " << out.start_pos 
             << ", 'length': " << out.length 
             << ", 'repeat_count': " << out.repeat_count 
             << ", 'is_inversion': " << (out.is_inversion ? "True" : "False") 
             << "}" << endl;
    }

    return 0;
}