//function for sorting a series of entries according to degrees in those entries
template<typename degree_type>
void sort_deg(std::vector<std::pair<matrix_index,matrix_index>> &row_cols, std::vector<degree_type> row_degs){
	//comparison function for entries according to the degree
	std::function<bool(const std::pair<matrix_index,matrix_index>&, const std::pair<matrix_index,matrix_index>&)> cmp = [&row_degs](const std::pair<matrix_index,matrix_index>& a, const std::pair<matrix_index,matrix_index>& b){ 
		return row_degs[a.first] < row_degs[b.first]; };
		//stable sorting the entries
		std::stable_sort(row_cols.begin(),row_cols.end(), cmp);
}
