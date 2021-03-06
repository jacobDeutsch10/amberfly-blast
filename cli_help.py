from PyInquirer import Validator, ValidationError
import os
class DirectoryValidator(Validator):
    def validate(self, doc):
        try:
            if not os.path.exists(str("./"+ doc.text)):
                raise ValidationError(
                    message='The directory at this location does not exist',
                    cursor_position=len(doc.text))  # Move cursor to end
        except ValueError:
            raise ValidationError(
                message='Please enter a valid directory path',
                cursor_position=len(doc.text))  # Move cursor to end

def get_numeric_input(msg):
	done = False
	while not done:
		num = input(msg)
		done = num.isnumeric()
		if not done:
			print("must be numeric")
	return int(num)

init_choices =  [
            'Do blast analysis from Patmos output',
            'Analyze previously created blast output', 
        ]
init_question = [
    {
        'type': 'list',
        'name': 'mode',
        'message': 'pick one',
        'choices': init_choices
    }
]
questions = [
    {
        'type': 'input',
        'name': 'PATMOS',
        'message': 'Enter location of PATMOS outputfile',
        'validate': DirectoryValidator,
        'default': 'singleFULL_INFO.xls'
    },
    {
        'type': 'input',
        'name': 'DB',
        'message': 'What would you like to name this database?',
    },
  
    {
        'type': 'checkbox',
        'qmark': '?',
        'message': 'Select behaviors',
        'name': 'behaviors',
        'choices': [ 
            {
                'name': 'dt'
            },
            {
                'name': 'xt'
            },
            {
                'name': 'yt'
            },
            {
                'name': 'zt',
            },
            {
                'name': 'xt_avg'
            },
            {
                'name': 'yt_avg'
            },
            {
                'name': 'zt_avg'
            },
            {
                'name': 'xxt'
            },
            {
                'name': 'xpt'
            },
            {
                'name': 'xp'
            },
            {
                'name': 'wd',
            },
            {
                'name': 'xp_avg'
            },
            {

                'name': 'wd_avg'
            },
            {
                'name': 'PCAsigma',
            },
            {
                'name': 'PCAmu'
            },
            {
                'name': 'PCAwt'
            },
            {
                'name': 'vxt'
            },
            {
                'name': 'vyt',
            },
            {
                'name': 'vzt'
            },
            {
                'name': 'vxt_avg'
            },
            {
                'name': 'vyt_avg'
            },
            {
                'name': 'vzt_avg'
            },
            {
                'name': 'vavg'
            },
            {
                'name': 'vv'
            },
            {
                'name': 'vv_avg'
            },
            {
                'name': 'vvt'
            },
            {
                'name': 'vpvy_avg'
            },
            {
                'name': 'axt'
            },
            {
                'name': 'ayt'
            },
            {
                'name': 'azt'
            },
            { 
                'name': 'aavg' 
            }, 
            {
                'name': 'axt_avg'
            },
            {
                'name': 'ayt_avg'
            },
            {
                'name': 'azt_avg'
            },
            {
                'name': 'curvature'
            }, 
            {
                'name': 'curveRadius'
            },
            {
                'name': 'vxp'
            }, 
            {
                'name': 'vxp_avg'
            },
            {
                'name': 'axp_47'
            },
            {
                'name': 'axp_avg'
            }, 
            { 
                'name': 'thetaToR'
            },
            {
                'name': 'psiToR'
            }, 
            {
                'name': 'aa'
            }
        ],
        'validate': lambda answer: 'You must choose at least one behavior' \
            if len(answer) == 0 else True
    }
]
