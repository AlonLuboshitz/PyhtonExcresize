# Alon Luboshitz 312115090
import sys
import pandas as pd
'''class data gets three paths for csv files - books,user,rating
validating three paths for that files holding memebers as data frames for those files.
validating spesific colums in those files for functions later to be done and keeping
those colums in a dict for each file.'''
class myData:
    def __init__(self,books_path,ratings_path,users_path) -> None:
        self.__dict_path = {}
        answer = self.valid_path([books_path,ratings_path,users_path]) 
        assert answer[0],"missing files: {}".format(answer[1])
        self.__books_data = pd.read_csv(self.__dict_path['book'],sep=";",encoding='latin-1',on_bad_lines='skip')
        self.__ratings_data = pd.read_csv(self.__dict_path['rating'],sep=";",encoding='latin-1',on_bad_lines='skip')
        self.__users_data = pd.read_csv(self.__dict_path['user'],sep=";",encoding='latin-1',on_bad_lines='skip')
        self.book_matching_dict = {}
        self.rating_matching_dict = {}
        self.user_matching_dict = {}
        valid_colums = self.validate_files()
        for bool in valid_colums:
            assert bool[0],'coudlnt find corresponding titles in the {} dataframe'.format(bool[1])
        self.remove_bad_lines() 
    '''this function check that the paths given to the ctor are for correspoding
    files. meaning there are three paths for books,user,rating.
    if not returning false and the missing files as a string'''
    def valid_path(self,path_list):
        files_text = ['book','user','rating']
        for text in files_text:
            for path in path_list:
                if text in path.lower():
                    self.__dict_path[text] = path
                    break
                else: self.__dict_path[text] = False
        missing_files = ''
        return_val = True
        for key, val in self.__dict_path.items():
            if not val:
                missing_files = missing_files + key + ','
                return_val = False
        return (return_val,missing_files)

    '''this function checks the file consits the disered data by checking the index line
    for books validate there is isbn,title,year,author.
    for rankings - isbn,id,rating
    for users - id, location,age
    returns a tuple of 3 with boolean expressions for books,ranking,users'''
    def validate_files(self):
        #check book
        books_colums_tocheck = ['isbn','title','year','author']
        books_colums = self.__books_data.columns.tolist()
        self.book_matching_dict = {text: [] for text in books_colums_tocheck}
        for text in books_colums_tocheck:
            for title in books_colums:
                if text in title.lower():
                    self.book_matching_dict[text] = title
                    break
        book_empty_list = (all(bool(lst) for lst in self.book_matching_dict.values()),'Books')
        #check rating
        rating_colums_tocheck = ['isbn','id','rating']
        rating_colums = self.__ratings_data.columns.tolist()
        self.rating_matching_dict = {text: [] for text in rating_colums_tocheck}
        for text in rating_colums_tocheck:
            for title in rating_colums:
                if text in title.lower():
                    self.rating_matching_dict[text] = title
                    break
        rating_empty_list = (all(bool(lst) for lst in self.rating_matching_dict.values()),'Rating')
        #check user
        user_colums_tocheck = ['location','id','age']
        user_colums = self.__users_data.columns.tolist()
        self.user_matching_dict = {text: [] for text in user_colums_tocheck}
        for text in user_colums_tocheck:
            for title in user_colums:
                if text in title.lower():
                    self.user_matching_dict[text] = title
                    break
        user_empty_list = (all(bool(lst) for lst in self.user_matching_dict.values()),'Users')
        return (book_empty_list,rating_empty_list,user_empty_list)
    '''this function remove lines with non numeric values in the year colum for the book df.
    it changes afterward the colum to an int type'''
    def remove_bad_lines(self):
        year = self.book_matching_dict['year']
        self.__books_data[year] = pd.to_numeric(self.__books_data[year],errors='coerce')
        self.__books_data = self.__books_data.dropna(subset=[year])
        self.__books_data.loc[:,year] = self.__books_data.loc[:,year].astype('int32')
    '''checks years are bigger then 0 and not decimal'''
    def check_years(self,year_list):
        for year in year_list:
            if not isinstance(year,int):
                return (False,'non integer value - {}'.format(year))
            if year <= 0:
                return (False,'bad year value - {} '.format(year))
        return (True,)
    def strip_colums(self):
       self.__ratings_data.loc[:,"ISBN"] = self.__ratings_data.loc[:,"ISBN"].str.strip().astype('str')
    def num_year(self,x,y):
        assert x != y, 'first year argument  - {} equals  second - {}'.format(x,y)
        assert x < y,('first year argument  - {} is bigger then second - {}'.format(x,y))
        answer = self.check_years([x,y]) 
        assert answer[0],answer[1]
        #filter true values of year range
        year = self.book_matching_dict['year']
        is_within_range = (self.__books_data[year] >= x) & (self.__books_data[year] < y)
        book_in_range = self.__books_data[is_within_range]
        return(len(book_in_range[year]))
    '''return spesific data frame with only books written in that year
    coloms returned are - title and author'''
    def df_published(self,year):
        answer = self.check_years([year]) 
        assert answer[0],answer[1]
        #filter only year values
        _year = self.book_matching_dict['year']
        in_year = (self.__books_data[_year] == year)
        author_title = pd.DataFrame(self.__books_data[in_year],
                                     columns=[self.book_matching_dict['title'],self.book_matching_dict['author'] ])
        return(author_title)
    '''return tuples list including number of books for each year
     sorted in ascending order '''
    #VALID X ND Y
    def num_books_by_year(self,x,y):
        answer = self.check_years([x,y]) 
        assert answer[0],answer[1]
        assert x <= y,('first year argument  - {} is bigger then second - {}'.format(x,y))
        #filter by years
        year = self.book_matching_dict['year']
        in_range_year_df = self.__books_data[(self.__books_data[year] >= x) & (self.__books_data[year] <= y)]
        count = in_range_year_df[year].astype('int32').value_counts().sort_index()
        tuples_list = list(count.items())
        return (tuples_list)
    '''given a country name function returns the mean age and std.
    cleaning null inputs in age colum.
    creating country colum by grouped countries with aggregation for std,mean build in pands funcitons
    try to locate the country, if there is no key returns string error'''
    def mean_std(self,country):
        location = self.user_matching_dict['location']
        age = self.user_matching_dict['age']
        #ADD TRY AND EXCEPT FOR COUNTRY
        # Replace NULL values with NaN and remove them
        self.__users_data[age] = self.__users_data[age].replace('NULL', pd.NA)
        self.__users_data = self.__users_data.dropna(subset=[age])
        #create new country colum with countries 

        #if cleaning the string further with those methods we get more countries!
        # .str.replace('.','').str.replace('"','')
        self.__users_data['Country'] = (self.__users_data.loc[:,location].str.split(',').str[-1].str.strip().astype('str'))
        #new DF grouped by countries with std, mean rounded to 3 digits.
        country_data = self.__users_data.groupby('Country').agg({age:['mean','std']})
        country_data = country_data.round(3)
        
        try:
            mean_std_tuple = tuple(country_data.loc[country,age])
            return mean_std_tuple
        except KeyError:
            return("No such country: {}, exsits in data".format(country))
    '''function gets a book name, retrive it isbn list and extract the mean rating
    from the isbn list from the rating df.'''
    def mean_rating(self,book_name):
        title = self.book_matching_dict['title']
        isbn = self.book_matching_dict['isbn']
        try:
            #group by title and get the group of the desired book
            books_title = self.__books_data.groupby(title).get_group(book_name)
        except KeyError:
            return("No book: {} have been found".format(book_name))
        #list the isbns of that book
        book_isbn = books_title[isbn].tolist()
        # filter the ratings data via that list
        filter_rating = self.__ratings_data[self.__ratings_data.loc[:,isbn].isin(book_isbn)]
        rating = self.rating_matching_dict['rating']
        #calculate the mean ratings after rating data been filtered.
        mean_reating = filter_rating.loc[:,rating].mean()
        return (mean_reating)
    '''function get k number of books rating to be presented in ascending ordrer
    thourh thier rating value'''
    def top_k(self,k):
        assert isinstance(k,int),"{} is non interger value".format(k)
        assert int(k) > 0 ,"number of users - {} cannot be 0 or negative".format(k)
        self.strip_colums()
        rating = self.rating_matching_dict['rating']
        author = self.book_matching_dict['author']
        title = self.book_matching_dict['title']
        rating_isbn = self.rating_matching_dict['isbn']
        #group the books by isbn and for each group calculate mean rating.
        grouped_isbns = self.__ratings_data.groupby(rating_isbn)[rating].mean() 
        #add the mean rating to the books data by the isbn value
        rated_books = self.__books_data.join(grouped_isbns,on=rating_isbn)
        #discard NA values in the rating colum
        rated_books = rated_books.dropna(subset=[rating])
        #sort values first by rating - descending, author - asc, title - asc
        rated_books =  rated_books.sort_values([rating,author,title],ascending=[False,True,True])
        #if k is bigger then df all df is being returned.
        return rated_books.loc[:,[title,author,rating]].head(int(k))
    '''this function group the users, sorts them by how many books they read'''
    def most_active(self,user):
        assert isinstance(user,int),"{} is non interger value".format(user)
        assert int(user) > 0 ,"user number - {} cannot be 0 or negative".format(user)
        id = self.rating_matching_dict['id']
        #clean id colum from tracing spaces
        self.__ratings_data.loc[:,id].astype('str').str.strip()
        #groupby id and calculate number of rows for each id group
        group_id= self.__ratings_data.groupby(id).size()
        #sort by ascending size value
        group_id = group_id.sort_values(ascending=False)
        try:
            answer = (group_id.iloc[int(user)-1])
            return answer
        except IndexError:
            return "user number given- {}, is bigger then amount of users!".format(user)
            



    