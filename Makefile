check: 
	python setup.py sdist 
	twine check dist/* 
	
upload: 
	twine upload dist/* 

merge: 
	git remote update 
	git merge upsteam/main 

clean: 
	rm -rf dist raser.egg-info  


