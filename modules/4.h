//addition
template<typename index, typename R>
vectors<index,R> ModuleOp<index,R>::add(const vectors<index,R> &x, const vectors<index,R> &y){   
	if(isZero(x)) return y;
	if(isZero(y)) return x;
	
	vectors<index,R> result;
	result.dataArray.reserve(x.size()+y.size());
	
	//the iterators to go through x and y
	auto itx = x.dataArray.begin();
	auto ity = y.dataArray.begin();
	
	while(true){
		//if reach the end of x, then copy the rest of y
		if(itx==x.dataArray.end()){
			if(ity==y.dataArray.end())
				//this is the exit of the loop
				break;
			else{
				result.push(*ity);	
				ity++;
			}
		}
		//if reach end of y, then copy the rest of x
		else if(ity==y.dataArray.end()){
			result.push(*itx);	
			itx++;
		}
		//if the current term in x is less, then copy it
		else if(itx->ind < ity->ind){
			result.push(*itx);	
			itx++;
		}
		//if the current term in y is less, then copy it
		else if(itx->ind > ity->ind){
			result.push(*ity);	
			ity++;
		}
		//if the current terms are equal, then add them
		else{
			typename vectors<index,R>::term newTerm(itx->ind, std::move(ringOper -> add(itx->coeficient, ity->coeficient)));
			if(!ringOper->isZero(newTerm.coeficient))
				result.push(std::move(newTerm));
			itx++;
			ity++;
		}
	}
	return result;
}

//destructive addition
template<typename index, typename R>
vectors<index,R> ModuleOp<index,R>::add(vectors<index,R> &&x, vectors<index,R>&&y){
	if(isZero(x)) return std::move(y);
	if(isZero(y)) return std::move(x);
	
	vectors<index,R> result;
	result.dataArray.reserve(x.size()+y.size());
	
	auto itx = x.dataArray.begin();
	auto ity = y.dataArray.begin();
	
	while(true){
		if(itx==x.dataArray.end()){
			if(ity==y.dataArray.end())
				break;
			else {
				result.push(std::move(*ity));
				ity++;
			}
		}
		else if(ity==y.dataArray.end()){
			result.push(std::move(*itx));	
			itx++;
		}
		else if(itx->ind < ity->ind){
			result.push(std::move(*itx));	
			itx++;
		}
		else if(itx->ind > ity->ind){
			result.push(std::move(*ity));	
			ity++;
		}
		else{
			typename vectors<index,R>::term newTerm(itx->ind, std::move(ringOper -> add(std::move(itx->coeficient), std::move(ity->coeficient))));
			if(!ringOper->isZero(newTerm.coeficient))
				result.push(std::move(newTerm));
			itx++;
			ity++;
		}
	}
	return result;
}
